"""Assemble the FS³ process (unit operations + connections) for a given recipe."""

import math
from bisect import bisect_right
from typing import Dict, List, Sequence, Tuple

import numpy as np
from numba import cfunc, float64

import fs3
from recipe import FRACTIONS, FlowSheetConfig, RecipeStep

# Pump percentage -> volumetric flow rate [m³/s] (mean of measured upper/lower calibration curves).
_PUMP_PERCENTAGES = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90], dtype=np.float64)
_UPPER_CURVE = np.array([0, 91, 220, 330, 465, 621, 783, 916, 1075, 1236, 1420, 1495, 1750, 1830, 1858, 1865, 1860, 1865, 1895], dtype=np.float64)
_LOWER_CURVE = np.array([0, 117, 260, 380, 505, 654, 783, 920, 1080, 1238, 1424, 1501, 1758, 1838, 1861, 1875, 1868, 1875, 1910], dtype=np.float64)
_FLOW_CURVE = (_UPPER_CURVE + _LOWER_CURVE) / 2.0 / 6e7

# Mixing percentage -> axial dispersion coefficient in the process chamber.
_MIXING_PERCENTAGES = np.array([0, 20, 40, 60, 80], dtype=np.float64)
_D_AX_CURVE = np.array([0.0, 1.3e-5, 5e-5, 6.9e-5, 8.4e-5], dtype=np.float64)


# Uncapture is ramped in over this many seconds at the start of a mixing step.
_UNCAPTURE_RAMP_S = 10.0
_EPS = 1e-9  # offset at which the left-hand limit of a step is sampled


@cfunc(float64(float64), cache=True)
def _a_eff_function(um_u0_ratio):
    a1, a2, a3, p, q = 2.035, 107.1, -0.00808, 0.07477, 1.083
    return 1.0 / (1.0 + a1 * math.exp(-p * um_u0_ratio) + a2 * math.exp(-q * um_u0_ratio) + a3)


def _table(f, boundaries: Sequence[float], extra_points: Sequence[float] = ()) -> fs3.LinearInterpolator:
    """Sample a function of time into a lookup table. Each point is stored twice -- left-hand
    limit and value -- so a repeated x jumps with zero width and steps stay exact at any time
    step. Pass ``extra_points`` where f is not constant within a section."""
    t = np.unique(np.concatenate([[0.0], boundaries, extra_points]))
    y = np.empty(2 * t.size)
    y[0::2] = [f(t_i - _EPS) for t_i in t]
    y[1::2] = [f(t_i) for t_i in t]
    return fs3.LinearInterpolator(np.repeat(t, 2)[1:], y[1:])


def build_process(
    rs, cs, recipe: Sequence[RecipeStep], solutions: Dict[str, np.ndarray], discretization_factor: float
) -> Tuple[fs3.Process, Dict[str, object]]:
    """Assemble the process. Returns ``(process, unit_operations)`` where a unit operation is a
    fraction iff its name is in ``FRACTIONS``."""
    water = solutions["Water (B5)"]

    cumulative_ends: List[float] = []
    total = 0.0
    for step in recipe:
        total += step.t_duration
        cumulative_ends.append(total)

    def interp_pump(p):
        return float(np.interp(p, _PUMP_PERCENTAGES, _FLOW_CURVE))

    def capture(t):
        """Capture flag (1.0) or uncapture rate; ramped in over the first seconds of a mixing step."""
        i = bisect_right(cumulative_ends, t)
        if i >= len(recipe) or recipe[i].mixing_percentage == 0.0:
            return 1.0
        rate = math.log(2.0) / 5.0
        ramp = (t - (cumulative_ends[i] - recipe[i].t_duration)) / _UNCAPTURE_RAMP_S
        return -rate * min(ramp * ramp, 1.0)

    # ----- Schedules the solver calls per right-hand side evaluation -----
    # Built as native lookups; a Python callable here would cost ~1 µs per call.
    section_starts = np.array([0.0, *cumulative_ends])

    def steps(value_of, tail=0.0):
        """Piecewise-constant schedule: value_of(step) within a step, tail after the recipe."""
        return fs3.PiecewiseConstantInterpolator(section_starts, np.array([*map(value_of, recipe), tail]))

    flow_rate = steps(lambda s: interp_pump(s.pump_percentage))
    pc_d_ax = steps(lambda s: float(np.interp(s.mixing_percentage, _MIXING_PERCENTAGES, _D_AX_CURVE)))
    flow_rate_line = steps(lambda s: interp_pump(s.pump_percentage) * (s.flow_sheet_configuration == FlowSheetConfig.LINE))
    flow_rate_loop = steps(lambda s: interp_pump(s.pump_percentage) * (s.flow_sheet_configuration == FlowSheetConfig.LOOP))
    flow_rate_fractions = {n: steps(lambda s, n=n: interp_pump(s.pump_percentage) * (s.fraction == n)) for n in FRACTIONS}
    inlet_concentration = fs3.PiecewiseConstantInterpolator(
        section_starts, np.array([solutions[step.inlet] for step in (*recipe, recipe[-1])])
    )

    # capture() varies within a step, so its ramp is sampled instead.
    ramp_points = [
        t
        for step, end in zip(recipe, cumulative_ends)
        if step.mixing_percentage
        for t in np.linspace(end - step.t_duration, end - step.t_duration + _UNCAPTURE_RAMP_S, 101)
    ]
    capture = _table(capture, cumulative_ends, ramp_points)

    # ----- Unit operations -----
    n_cells = lambda base: max(1, int(round(base * discretization_factor)))
    d_pipes = 0.0096
    a_cross_pipes = math.pi / 4.0 * d_pipes * d_pipes
    l_inlet_pipe = 0.67 + (0.00005 / a_cross_pipes)

    inlet = fs3.Inlet(rs, inlet_concentration)
    pipe_inlet = fs3.Pipe(rs, n_cells(10), a_cross_pipes, l_inlet_pipe, flow_rate, 0.0)
    pipe_inlet.set_const_initial_concentration(water)

    d_pc, d_shaft_pc, l_pc = 9.031e-2, 6.691e-2, 9.755e-02
    a_cross_pc = (math.pi / 4.0) * (d_pc * d_pc - d_shaft_pc * d_shaft_pc)
    porosity_pc = 2.5558e-4 / (a_cross_pc * l_pc)
    pc = fs3.RS_MagneticCaptureProcessChamber(
        rs, n_cells(30), a_cross_pc, l_pc, porosity_pc, 25, 8.0e-4, 19, 2,
        40e-3, 2.227e5, _a_eff_function, capture, flow_rate, pc_d_ax,
    )
    pc.set_const_initial_concentration(water)

    dead_volume = fs3.Volume(rs, 4.113e-5 + 8.381e-5)
    dead_volume.set_const_initial_concentration(water)

    pipe_outlet = fs3.Pipe(rs, n_cells(10), a_cross_pipes, 1.02, flow_rate, 0.0)
    pipe_outlet.set_const_initial_concentration(water)

    pipe_loop = fs3.Pipe(rs, n_cells(10), a_cross_pipes, 1.73, flow_rate_loop, 0.0)
    pipe_loop.set_const_initial_concentration(water)

    fractions = {name: fs3.Volume(rs) for name in FRACTIONS}

    # ----- Connections -----
    process = fs3.Process(cs, [inlet, pipe_inlet, pc, dead_volume, pipe_outlet, pipe_loop, *fractions.values()])
    process.add_connection(inlet.exit(), pipe_inlet.entry(), flow_rate_line)
    process.add_connection(pipe_inlet.exit(), pc.entry(), flow_rate)
    process.add_connection(pc.exit(), dead_volume.entry(), flow_rate)
    process.add_connection(dead_volume.exit(), pipe_outlet.entry(), flow_rate)
    for name, frac in fractions.items():
        process.add_connection(pipe_outlet.exit(), frac.entry(), flow_rate_fractions[name])
    process.add_connection(pipe_outlet.exit(), pipe_loop.entry(), flow_rate_loop)
    process.add_connection(pipe_loop.exit(), pipe_inlet.entry(), flow_rate_loop)

    unit_operations = {
        "inlet": inlet,
        "pipe_inlet": pipe_inlet,
        "pc": pc,
        "dead_volume": dead_volume,
        "pipe_outlet": pipe_outlet,
        "pipe_loop": pipe_loop,
        **fractions,
    }
    return process, unit_operations

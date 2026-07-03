"""Assemble the FS³ process (unit operations + connections) for a given recipe."""

import math
from bisect import bisect_right
from typing import Dict, List, Sequence, Tuple

import numpy as np

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


def _a_eff_function(um_u0_ratio: float) -> float:
    a1, a2, a3, p, q = 2.035, 107.1, -0.00808, 0.07477, 1.083
    return 1.0 / (1.0 + a1 * math.exp(-p * um_u0_ratio) + a2 * math.exp(-q * um_u0_ratio) + a3)


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

    def section_idx(t):
        i = bisect_right(cumulative_ends, t)
        return None if i >= len(recipe) else i

    def interp_pump(p):
        return float(np.interp(p, _PUMP_PERCENTAGES, _FLOW_CURVE))

    def flow_rate(t):
        i = section_idx(t)
        return interp_pump(0.0 if i is None else recipe[i].pump_percentage)

    def pc_d_ax(t):
        i = section_idx(t)
        return float(np.interp(0.0 if i is None else recipe[i].mixing_percentage, _MIXING_PERCENTAGES, _D_AX_CURVE))

    def capture(t):
        i = section_idx(t)
        if i is None:
            return 1.0
        step = recipe[i]
        if step.mixing_percentage == 0.0:
            return 1.0
        t_in_section = t - (cumulative_ends[i] - step.t_duration)
        uncapture_rate = math.log(2.0) / 5.0
        if t_in_section < 10.0:
            ramp = t_in_section / 10.0
            return -(ramp * ramp) * uncapture_rate
        return -uncapture_rate

    def inlet_concentration(t):
        i = section_idx(t)
        return solutions[(recipe[-1] if i is None else recipe[i]).inlet]

    def flow_rate_for(t, config):
        i = section_idx(t)
        if i is None:
            return interp_pump(0.0)
        return interp_pump(recipe[i].pump_percentage) if recipe[i].flow_sheet_configuration == config else 0.0

    def flow_rate_line(t):
        return flow_rate_for(t, FlowSheetConfig.LINE)

    def flow_rate_loop(t):
        return flow_rate_for(t, FlowSheetConfig.LOOP)

    def flow_rate_fraction(t, fraction_name):
        i = section_idx(t)
        if i is None:
            return interp_pump(0.0)
        return interp_pump(recipe[i].pump_percentage) if recipe[i].fraction == fraction_name else 0.0

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
        process.add_connection(pipe_outlet.exit(), frac.entry(), lambda t, n=name: flow_rate_fraction(t, n))
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

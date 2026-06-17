"""Python port of src/HGMS_hIgG_separation.cpp using FS³'s python bindings."""

import math
from bisect import bisect_right
from datetime import datetime
from pathlib import Path
from typing import Optional, Sequence, Tuple, Union

import numpy as np

import fs3
from model import build_component_system, build_reaction_system, build_solutions
from shared_data import FRACTIONS, FlowSheetConfig, RecipeStep, default_recipe


def run_HGMS_hIgG_separation(
    kf_ion: float,
    tau_reaction: float,
    save_obs: bool,
    timeout_seconds: float,
    solver_type,
    discretization_factor: float,
    recipe_steps: Optional[Sequence[RecipeStep]] = None,
    output_base_dir: Optional[Union[str, Path]] = None,
    run_tag: Optional[str] = None,
    return_paths: bool = False,
):
    cs = build_component_system()
    rs, activity_model = build_reaction_system(cs, tau_reaction)
    solutions = build_solutions(cs, rs, solver_type, timeout_seconds)
    water = solutions["Water (B5)"].concentrations

    recipe = list(recipe_steps) if recipe_steps is not None else default_recipe(solutions)

    # =============================================================
    # ========================== RECIPE ===========================
    # =============================================================
    fractions = {name: fs3.Volume(rs) for name in FRACTIONS}

    pump_percentages = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90], dtype=np.float64)
    upper_curve = np.array([0, 91, 220, 330, 465, 621, 783, 916, 1075, 1236, 1420, 1495, 1750, 1830, 1858, 1865, 1860, 1865, 1895], dtype=np.float64)
    lower_curve = np.array([0, 117, 260, 380, 505, 654, 783, 920, 1080, 1238, 1424, 1501, 1758, 1838, 1861, 1875, 1868, 1875, 1910], dtype=np.float64)
    flow_curve = (upper_curve + lower_curve) / 2.0 / 6e7

    mixing_percentages = np.array([0, 20, 40, 60, 80], dtype=np.float64)
    d_ax_curve = np.array([0.0, 1.3e-5, 5e-5, 6.9e-5, 8.4e-5], dtype=np.float64)

    cumulative_ends = []
    total_duration = 0.0
    for step in recipe:
        total_duration += step.t_duration
        cumulative_ends.append(total_duration)

    def section_idx(t: float):
        i = bisect_right(cumulative_ends, t)
        return None if i >= len(recipe) else i

    def interp_pump(p: float) -> float:
        return float(np.interp(p, pump_percentages, flow_curve))

    def interp_dax(p: float) -> float:
        return float(np.interp(p, mixing_percentages, d_ax_curve))

    def flowRateFunction(t: float) -> float:
        i = section_idx(t)
        return interp_pump(0.0 if i is None else recipe[i].pump_percentage)

    def pc_D_ax_function(t: float) -> float:
        i = section_idx(t)
        return interp_dax(0.0 if i is None else recipe[i].mixing_percentage)

    def capture_function(t: float) -> float:
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

    def inlet_function(t: float):
        i = section_idx(t)
        return (recipe[-1] if i is None else recipe[i]).inlet.concentrations

    def flowRate_for_config(t: float, *configs) -> float:
        i = section_idx(t)
        if i is None:
            return interp_pump(0.0)
        return interp_pump(recipe[i].pump_percentage) if recipe[i].flow_sheet_configuration in configs else 0.0

    def flowRate_function_line(t: float) -> float:
        return flowRate_for_config(t, FlowSheetConfig.LINE)

    def flowRate_function_loop(t: float) -> float:
        return flowRate_for_config(t, FlowSheetConfig.LOOP)

    def flowRate_function_frac(t: float, fraction_name: str) -> float:
        i = section_idx(t)
        if i is None:
            return interp_pump(0.0)
        return interp_pump(recipe[i].pump_percentage) if recipe[i].fraction == fraction_name else 0.0

    # =============================================================
    # ========================== UNIT OPERATIONS ==================
    # =============================================================
    inlet = fs3.Inlet(rs, inlet_function)

    d_pipes = 0.0096
    a_cross_pipes = math.pi / 4.0 * d_pipes * d_pipes
    l_inlet_pipe = 0.67 + (0.00005 / a_cross_pipes)
    scale_n_cells = lambda base: max(1, int(round(base * discretization_factor)))

    pipe_inlet = fs3.Pipe(rs, scale_n_cells(10), a_cross_pipes, l_inlet_pipe, flowRateFunction, 0.0)
    pipe_inlet.set_const_initial_concentration(water)

    def a_eff_function(um_u0_ratio: float) -> float:
        a1 = 2.035
        a2 = 107.1
        a3 = -0.00808
        p = 0.07477
        q = 1.083
        return 1.0 / (1.0 + a1 * math.exp(-p * um_u0_ratio) + a2 * math.exp(-q * um_u0_ratio) + a3)

    l_pc = 9.755e-02
    d_pc = 9.031e-2
    d_shaft_pc = 6.691e-2
    a_cross_pc = (math.pi / 4.0) * (d_pc * d_pc - d_shaft_pc * d_shaft_pc)
    v_int_pc = 2.5558e-4
    porosity_pc = v_int_pc / (a_cross_pc * l_pc)
    magnetic_field_strength = 2.227e5
    mnp_capacity = 40e-3

    pc = fs3.RS_MagneticCaptureProcessChamber(
        rs,
        scale_n_cells(30),
        a_cross_pc,
        l_pc,
        porosity_pc,
        25,
        8.0e-4,
        19,
        2,
        mnp_capacity,
        magnetic_field_strength,
        a_eff_function,
        capture_function,
        flowRateFunction,
        pc_D_ax_function,
    )
    pc.set_const_initial_concentration(water)

    v_dead = 4.113e-5 + 8.381e-5
    dead_volume = fs3.Volume(rs, v_dead)
    dead_volume.set_const_initial_concentration(water)

    pipe_outlet = fs3.Pipe(rs, scale_n_cells(10), a_cross_pipes, 1.02, flowRateFunction, 0.0)
    pipe_outlet.set_const_initial_concentration(water)

    pipe_loop = fs3.Pipe(rs, scale_n_cells(10), a_cross_pipes, 1.73, flowRate_function_loop, 0.0)
    pipe_loop.set_const_initial_concentration(water)

    print(f"Total duration of the process: {total_duration} seconds.")

    process = fs3.Process(cs, [inlet, pipe_inlet, pc, dead_volume, pipe_outlet, pipe_loop, *fractions.values()])

    process.add_connection(inlet.exit(), pipe_inlet.entry(), flowRate_function_line)
    process.add_connection(pipe_inlet.exit(), pc.entry(), flowRateFunction)
    process.add_connection(pc.exit(), dead_volume.entry(), flowRateFunction)
    process.add_connection(dead_volume.exit(), pipe_outlet.entry(), flowRateFunction)

    for name, frac in fractions.items():
        process.add_connection(pipe_outlet.exit(), frac.entry(), lambda t, n=name: flowRate_function_frac(t, n))

    process.add_connection(pipe_outlet.exit(), pipe_loop.entry(), flowRate_function_loop)
    process.add_connection(pipe_loop.exit(), pipe_inlet.entry(), flowRate_function_loop)

    # =============================================================
    # ========================== OBSERVER =========================
    # =============================================================
    solver = fs3.Solver(process, solver_type)

    dt_obs = 2.0
    n_obs_time_steps = int(total_duration / dt_obs) + 1
    print(f"Number of observation time steps: {n_obs_time_steps}")

    observed = {
        "pipe_inlet": pipe_inlet.all(),
        "pc_liquid": pc.liquid(),
        "pc_slurry": pc.slurry(),
        "pipe_outlet": pipe_outlet.all(),
        "pipe_loop": pipe_loop.all(),
        **{name: frac.all() for name, frac in fractions.items()},
    }
    observers = {
        name: fs3.TimeSeriesObserver(0.0, total_duration, n_obs_time_steps, target, solver, save_obs, True)
        for name, target in observed.items()
    }
    obs_pc_outlet = fs3.TimeSeriesObserver(0.0, total_duration, n_obs_time_steps, pc.exit(), solver, save_obs)

    # =============================================================
    # ========================== SOLVING ==========================
    # =============================================================
    solver.solve(total_duration, timeout_seconds)
    t_solve_duration = solver.get_solve_time()

    # =============================================================
    # ========================== RESULTS ==========================
    # =============================================================
    internal_time_stamps = solver.get_internal_time_stamps()

    output_base_dir = Path.cwd() if output_base_dir is None else Path(output_base_dir)
    run_name = run_tag or datetime.now().strftime("run_%Y-%m-%d_%H-%M-%S")
    run_dir = output_base_dir / run_name
    obs_dir = run_dir / "obs"
    if save_obs:
        (run_dir / "logs").mkdir(parents=True, exist_ok=True)
        (run_dir / "bench").mkdir(parents=True, exist_ok=True)
        obs_dir.mkdir(parents=True, exist_ok=True)

        np.save(obs_dir / "internal_timestamps.npy", np.asarray(internal_time_stamps, dtype=np.float64))

        obs_unit_ops = str(obs_dir / "unit_operations.npz")
        for name, obs in observers.items():
            obs.save_to_npz(obs_unit_ops, name, "a")

        obs_pc_outlet.save_to_npy(str(obs_dir / "pc_outlet.npy"))

        cell_index = pipe_outlet.n_cells // 2
        pH_activity = fs3.convert_to_pH(observers["pipe_outlet"], cell_index, pipe_outlet.cell_volume, cs, activity_model)
        np.save(obs_dir / "pipe_outlet_middle_cell_pH_activity.npy", np.asarray(pH_activity, dtype=np.float64))

        pH_concentration = fs3.convert_to_pH(observers["pipe_outlet"], cell_index, pipe_outlet.cell_volume, cs, fs3.NoActivityModel(cs))
        np.save(obs_dir / "pipe_outlet_middle_cell_pH_concentration.npy", np.asarray(pH_concentration, dtype=np.float64))

    print(f"Solved with {len(internal_time_stamps)} internal time stamps in {t_solve_duration} seconds.")
    if t_solve_duration > timeout_seconds:
        print(f"WARNING: Solver timed out after {timeout_seconds} seconds at t={solver.get_t()} / {total_duration}")

    result: Tuple[float, float, int] = (float("nan"), t_solve_duration, len(internal_time_stamps))
    if return_paths:
        return result, run_dir, obs_dir
    return result


if __name__ == "__main__":
    run_HGMS_hIgG_separation(
        1e3,
        0.1,
        True,
        float("inf"),
        fs3.SolverType.ERK,
        0.5,
    )

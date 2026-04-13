"""Single-function Python port of src/HGMS_hIgG_separation.cpp using FS3 bindings.<
Currently it has no lambda functions and therefore no Gouy-Chapman model and Langmuir absorption."""

from __future__ import annotations

import math
from bisect import bisect_right
from datetime import datetime
from pathlib import Path

import numpy as np
import fs3


def run_HGMS_hIgG_separation(
    kf_ion: float,
    tau_reaction: float,
    save_obs: bool,
    timeout_seconds: float,
    solver_type,
    discretization_factor: float,
):
    del kf_ion  # Kept for signature parity with C++.

    NO_FLOW = 0
    LINE = 1
    LOOP = 2

    # =============================================================
    # ========================== COMPONENTS =======================
    # =============================================================
    cs = fs3.ComponentSystem(
        [
            fs3.Component("H₂O", molar_mass_kg_per_mol=18.01528e-3),
            fs3.Component("H⁺", charge=1, molar_mass_kg_per_mol=1.008e-3, truesdell_jones_alpha=4.78e-10, truesdell_jones_beta=0.24e-3),
            fs3.Component("OH⁻", charge=-1, molar_mass_kg_per_mol=17.007e-3, truesdell_jones_alpha=10.65e-10, truesdell_jones_beta=0.21e-3),
            fs3.Component("Na⁺", charge=1, molar_mass_kg_per_mol=22.990e-3, truesdell_jones_alpha=4.32e-10, truesdell_jones_beta=0.06e-3),
            fs3.Component("Cl⁻", charge=-1, molar_mass_kg_per_mol=35.453e-3, truesdell_jones_alpha=3.71e-10, truesdell_jones_beta=0.01e-3),
            fs3.Component("TrisH⁺", charge=1, molar_mass_kg_per_mol=121.14e-3, truesdell_jones_alpha=4.0e-10, truesdell_jones_beta=0.0e-3),
            fs3.Component("Tris", molar_mass_kg_per_mol=121.14e-3),
            fs3.Component("AcH", molar_mass_kg_per_mol=60.052e-3),
            fs3.Component("Ac⁻", charge=-1, molar_mass_kg_per_mol=59.044e-3, truesdell_jones_alpha=4.5e-10, truesdell_jones_beta=0.0e-3),
            fs3.Component("GlyH₂⁺", charge=1, molar_mass_kg_per_mol=75.067e-3),
            fs3.Component("GlyH", molar_mass_kg_per_mol=75.067e-3),
            fs3.Component("Gly⁻", charge=-1, molar_mass_kg_per_mol=74.059e-3, truesdell_jones_alpha=4.0e-10, truesdell_jones_beta=0.0e-3),
            fs3.Component("hIgG", molar_mass_kg_per_mol=150.0),
            fs3.Component("MNP-OH₂⁺", type=fs3.ComponentType.MagneticNanoParticleGroup, molar_mass_kg_per_mol=18.015e-3),
            fs3.Component("MNP-OH", type=fs3.ComponentType.MagneticNanoParticleGroup, molar_mass_kg_per_mol=17.007e-3),
            fs3.Component("MNP-O⁻", type=fs3.ComponentType.MagneticNanoParticleGroup, molar_mass_kg_per_mol=15.999e-3),
            fs3.Component("MNP-hIgG", type=fs3.ComponentType.MagneticNanoParticleGroup, molar_mass_kg_per_mol=150000e-3),
            fs3.Component("MNP1", type=fs3.ComponentType.MagneticNanoParticle, radius_m=1039e-9, density_kg_per_m3=5170, magnetic_saturation_A_per_m=3.5e5),
            fs3.Component("MNP2", type=fs3.ComponentType.MagneticNanoParticle, radius_m=1209e-9, density_kg_per_m3=5170, magnetic_saturation_A_per_m=3.5e5),
            fs3.Component("MNP3", type=fs3.ComponentType.MagneticNanoParticle, radius_m=1406e-9, density_kg_per_m3=5170, magnetic_saturation_A_per_m=3.5e5),
            fs3.Component("MNP4", type=fs3.ComponentType.MagneticNanoParticle, radius_m=1635e-9, density_kg_per_m3=5170, magnetic_saturation_A_per_m=3.5e5),
            fs3.Component("MNP5", type=fs3.ComponentType.MagneticNanoParticle, radius_m=1901e-9, density_kg_per_m3=5170, magnetic_saturation_A_per_m=3.5e5),
            fs3.Component("MNP6", type=fs3.ComponentType.MagneticNanoParticle, radius_m=2211e-9, density_kg_per_m3=5170, magnetic_saturation_A_per_m=3.5e5),
            fs3.Component("MNP7", type=fs3.ComponentType.MagneticNanoParticle, radius_m=2571e-9, density_kg_per_m3=5170, magnetic_saturation_A_per_m=3.5e5),
            fs3.Component("MNP8", type=fs3.ComponentType.MagneticNanoParticle, radius_m=2990e-9, density_kg_per_m3=5170, magnetic_saturation_A_per_m=3.5e5),
            fs3.Component("MNP9", type=fs3.ComponentType.MagneticNanoParticle, radius_m=3477e-9, density_kg_per_m3=5170, magnetic_saturation_A_per_m=3.5e5),
            fs3.Component("MNP10", type=fs3.ComponentType.MagneticNanoParticle, radius_m=4043e-9, density_kg_per_m3=5170, magnetic_saturation_A_per_m=3.5e5),
        ]
    )

    # =============================================================
    # ========================== REACTIONS ========================
    # =============================================================
    activity_model = fs3.TruesdellJonesActivityModel(cs)
    rs = fs3.ReactionSystem(cs, activity_model)

    k_w = 10.0 ** (-14) * 1e6
    c_w = 1000.0 / cs["H₂O"].molar_mass

    ks_w = k_w / c_w
    ks_tris_h_plus = 10.0 ** (-8.3) * 1e3
    ks_ac_h = 10.0 ** (-4.75) * 1e3
    ks_gly_h2_plus = 10.0 ** (-2.35) * 1e3
    ks_gly_h = 10.0 ** (-9.78) * 1e3
    ks_mnp_oh2_plus = 10.0 ** (-(7.9 - 0.5 * 2.5)) * 1e3
    ks_mnp_oh = 10.0 ** (-(7.9 + 0.5 * 2.5)) * 1e3

    rs.add(fs3.mass_action_law_inverse_rate_prediction(cs, "H₂O <=> H⁺ + OH⁻", ks_w, tau_reaction, cs.get_idx("H⁺")))
    rs.add(fs3.mass_action_law_inverse_rate_prediction(cs, "TrisH⁺ <=> Tris + H⁺", ks_tris_h_plus, tau_reaction, cs.get_idx("H⁺")))
    rs.add(fs3.mass_action_law_inverse_rate_prediction(cs, "AcH <=> Ac⁻ + H⁺", ks_ac_h, tau_reaction, cs.get_idx("H⁺")))
    rs.add(fs3.mass_action_law_inverse_rate_prediction(cs, "GlyH₂⁺ <=> GlyH + H⁺", ks_gly_h2_plus, tau_reaction, cs.get_idx("H⁺")))
    rs.add(fs3.mass_action_law_inverse_rate_prediction(cs, "GlyH <=> Gly⁻ + H⁺", ks_gly_h, tau_reaction, cs.get_idx("H⁺")))

    # Closest available binding. The custom C++ Gouy-Chapman callback is not exposed in Python yet.
    rs.add(fs3.mass_action_law_inverse_rate_prediction(cs, "MNP-OH₂⁺ <=> MNP-OH + H⁺", ks_mnp_oh2_plus, tau_reaction, cs.get_idx("H⁺")))
    rs.add(fs3.mass_action_law_inverse_rate_prediction(cs, "MNP-OH <=> MNP-O⁻ + H⁺ ", ks_mnp_oh, tau_reaction, cs.get_idx("H⁺")))

    # =============================================================
    # ========================== SOLUTIONS ========================
    # =============================================================
    buffer_eps = 0.0
    t_reaction_duration = 100.0
    mnp_ns = 1.08e-5
    mnp_specific_surface = 1e5

    def one_cell_reaction(solution: np.ndarray, duration: float) -> tuple[np.ndarray, float]:
        volume = fs3.Volume(rs, 1.0)
        volume.set_const_initial_concentration(solution)
        process_local = fs3.Process(cs, [volume])
        solver_local = fs3.Solver(process_local, solver_type)
        observer = fs3.SnapshotObserver(duration, volume.all(), True)
        solver_local.add(observer)
        solver_local.solve(t_stop=duration, timeout_seconds=timeout_seconds)
        return np.asarray(observer.snapshot, dtype=np.float64).reshape(-1), float(observer.error)

    c_mnp_feed = 16.5
    v_mnps_feed = 0.000588
    c_higg_feed = 0.45
    v_higg_feed = 0.0046
    c_hcl_feed = 0.01835e3
    c_tris_feed = 0.02e3
    c_nacl_feed = 0.15e3
    v_buffer_feed = 0.000135 + 0.006 - 0.6e-3 - 0.873e-3
    v_feed = v_mnps_feed + v_higg_feed + v_buffer_feed

    c_mnp_feed *= v_mnps_feed / v_feed
    c_higg_feed *= v_higg_feed / v_feed
    c_hcl_feed *= v_buffer_feed / v_feed
    c_tris_feed *= v_buffer_feed / v_feed
    c_nacl_feed *= v_buffer_feed / v_feed

    feed = np.full(cs.n_components, buffer_eps, dtype=np.float64)
    feed[cs.get_idx("H₂O")] = c_w
    feed[cs.get_idx("H⁺")] = c_hcl_feed
    feed[cs.get_idx("Na⁺")] = c_nacl_feed
    feed[cs.get_idx("Cl⁻")] = c_hcl_feed + c_nacl_feed
    feed[cs.get_idx("Tris")] = c_tris_feed
    feed[cs.get_idx("hIgG")] = c_higg_feed

    mnp_distribution = np.array([0.04852, 0.1853, 0.25496, 0.14803, 0.03495, 0.03743, 0.09137, 0.10777, 0.068, 0.02367], dtype=np.float64)
    if abs(float(mnp_distribution.sum()) - 1.0) > 1e-3:
        raise RuntimeError("MNP distribution does not sum to 1.")
    for i in range(10):
        feed[cs.get_idx(f"MNP{i + 1}")] = c_mnp_feed * mnp_distribution[i]
    feed[cs.get_idx("MNP-OH")] = mnp_ns * mnp_specific_surface * c_mnp_feed

    feed, feed_error = one_cell_reaction(feed, t_reaction_duration * 10.0)
    pH_feed = -math.log10(feed[cs.get_idx("H⁺")] * 1e-3)
    q_feed = feed[cs.get_idx("MNP-hIgG")] / feed[cs.get_idx("MNP1") : cs.get_idx("MNP10") + 1].sum()
    q_max = float(np.interp(pH_feed, np.array([2.0, 2.5, 2.6, 3.5, 3.7, 4.5, 7.0]), np.array([0.00, 0.01, 0.02, 0.05, 0.13, 0.28, 0.31])))
    print(f"pH feed: {pH_feed} (should be ~7.3)")
    print(f"q feed: {q_feed} q_max at this pH: {q_max}")
    print(f"Feed reaction error: {feed_error}")

    buffer1 = np.full(cs.n_components, buffer_eps, dtype=np.float64)
    buffer1[cs.get_idx("H₂O")] = c_w
    buffer1[cs.get_idx("H⁺")] = 0.019156e3
    buffer1[cs.get_idx("Na⁺")] = 0.15e3
    buffer1[cs.get_idx("Cl⁻")] = 0.019156e3 + 0.15e3
    buffer1[cs.get_idx("Tris")] = 0.02e3
    buffer1, buffer1_error = one_cell_reaction(buffer1, t_reaction_duration)
    print(f"pH B1: {-math.log10(buffer1[cs.get_idx('H⁺')] * 1e-3)} (should be ~7)")
    print(f"B1 reaction error: {buffer1_error}")

    buffer2 = np.full(cs.n_components, buffer_eps, dtype=np.float64)
    buffer2[cs.get_idx("H₂O")] = c_w
    buffer2[cs.get_idx("H⁺")] = 0.072107e3
    buffer2[cs.get_idx("Na⁺")] = 0.2e3
    buffer2[cs.get_idx("Cl⁻")] = 0.072107e3
    buffer2[cs.get_idx("Ac⁻")] = 0.2e3
    buffer2, buffer2_error = one_cell_reaction(buffer2, t_reaction_duration)
    print(f"pH B2: {-math.log10(buffer2[cs.get_idx('H⁺')] * 1e-3)} (should be ~4.75)")
    print(f"B2 reaction error: {buffer2_error}")

    buffer3 = np.full(cs.n_components, buffer_eps, dtype=np.float64)
    buffer3[cs.get_idx("H₂O")] = c_w
    buffer3[cs.get_idx("H⁺")] = 0.20221e3
    buffer3[cs.get_idx("Na⁺")] = 0.2e3
    buffer3[cs.get_idx("Cl⁻")] = 0.20221e3
    buffer3[cs.get_idx("Ac⁻")] = 0.2e3
    buffer3, buffer3_error = one_cell_reaction(buffer3, t_reaction_duration)
    print(f"pH B3: {-math.log10(buffer3[cs.get_idx('H⁺')] * 1e-3)} (should be ~2.5)")
    print(f"B3 reaction error: {buffer3_error}")

    buffer5_water = np.full(cs.n_components, buffer_eps, dtype=np.float64)
    buffer5_water[cs.get_idx("H₂O")] = c_w
    buffer5_water, buffer5_water_error = one_cell_reaction(buffer5_water, t_reaction_duration)
    print(f"pH B5: {-math.log10(buffer5_water[cs.get_idx('H⁺')] * 1e-3)} (should be ~7.0 / 6.5)")
    print(f"B5 reaction error: {buffer5_water_error}")

    # =============================================================
    # ========================== RECIPE ===========================
    # =============================================================
    frac_feed_1 = fs3.Volume(rs)
    frac_feed_2 = fs3.Volume(rs)
    frac_wash_1 = fs3.Volume(rs)
    frac_wash_2 = fs3.Volume(rs)
    frac_wash_3 = fs3.Volume(rs)
    frac_elution_1 = fs3.Volume(rs)
    frac_elution_2 = fs3.Volume(rs)
    frac_elution_3 = fs3.Volume(rs)
    frac_elution_4 = fs3.Volume(rs)
    frac_elution_5 = fs3.Volume(rs)

    recipe = [
        ("Load 1", 1199, 20, 0, LINE, feed, frac_feed_1),
        ("Load 1 pause", 34, 0, 0, NO_FLOW, buffer5_water, None),
        ("Load 2", 96, 20, 0, LINE, feed, frac_feed_1),
        ("Load 2 pause", 33, 0, 0, NO_FLOW, buffer5_water, None),
        ("Wash 1 fill", 50, 40, 0, LINE, buffer1, frac_feed_2),
        ("Wash 1 resuspend", 65, 0, 60, NO_FLOW, buffer5_water, None),
        ("Wash 1 resuspend loop", 51, 40, 50, LOOP, buffer5_water, None),
        ("Wash 1 recapture", 69, 0, 0, NO_FLOW, buffer5_water, None),
        ("Wash 1 recapture loop", 196, 30, 0, LOOP, buffer5_water, None),
        ("Wash 1 pause", 6, 0, 0, NO_FLOW, buffer5_water, None),
        ("Wash 2 fill", 50, 40, 0, LINE, buffer1, frac_wash_1),
        ("Wash 2 resuspend", 65, 0, 60, NO_FLOW, buffer5_water, None),
        ("Wash 2 resuspend loop", 51, 40, 50, LOOP, buffer5_water, None),
        ("Wash 2 recapture", 69, 0, 0, NO_FLOW, buffer5_water, None),
        ("Wash 2 recapture loop", 196, 30, 0, LOOP, buffer5_water, None),
        ("Wash 2 pause", 14, 0, 0, NO_FLOW, buffer5_water, None),
        ("Wash 3 fill", 50, 40, 0, LINE, buffer2, frac_wash_2),
        ("Wash 3 resuspend", 65, 0, 60, NO_FLOW, buffer5_water, None),
        ("Wash 3 resuspend loop", 51, 40, 50, LOOP, buffer5_water, None),
        ("Wash 3 recapture", 69, 0, 0, NO_FLOW, buffer5_water, None),
        ("Wash 3 recapture loop", 196, 30, 0, LOOP, buffer5_water, None),
        ("Wash 3 pause", 13, 0, 0, NO_FLOW, buffer5_water, None),
        ("Elution 1 fill", 52, 40, 0, LINE, buffer3, frac_wash_3),
        ("Elution 1 resuspend", 90, 0, 60, NO_FLOW, buffer5_water, None),
        ("Elution 1 resuspend loop", 300, 30, 40, LOOP, buffer5_water, None),
        ("Elution 1 recapture", 60, 0, 0, NO_FLOW, buffer5_water, None),
        ("Elution 1 recapture loop", 200, 30, 0, LOOP, buffer5_water, None),
        ("Elution 1 pause", 7, 0, 0, NO_FLOW, buffer5_water, None),
        ("Elution 2 fill", 52, 40, 0, LINE, buffer3, frac_elution_1),
        ("Elution 2 resuspend", 90, 0, 60, NO_FLOW, buffer5_water, None),
        ("Elution 2 resuspend loop", 300, 30, 40, LOOP, buffer5_water, None),
        ("Elution 2 recapture", 60, 0, 0, NO_FLOW, buffer5_water, None),
        ("Elution 2 recapture loop", 200, 30, 0, LOOP, buffer5_water, None),
        ("Elution 2 pause", 9, 0, 0, NO_FLOW, buffer5_water, None),
        ("Elution 3 fill", 52, 40, 0, LINE, buffer3, frac_elution_2),
        ("Elution 3 resuspend", 90, 0, 60, NO_FLOW, buffer5_water, None),
        ("Elution 3 resuspend loop", 300, 30, 40, LOOP, buffer5_water, None),
        ("Elution 3 recapture", 60, 0, 0, NO_FLOW, buffer5_water, None),
        ("Elution 3 recapture loop", 200, 30, 0, LOOP, buffer5_water, None),
        ("Elution 3 pause", 10, 0, 0, NO_FLOW, buffer5_water, None),
        ("Elution 4 fill", 52, 40, 0, LINE, buffer3, frac_elution_3),
        ("Elution 4 resuspend", 90, 0, 60, NO_FLOW, buffer5_water, None),
        ("Elution 4 resuspend loop", 300, 30, 40, LOOP, buffer5_water, None),
        ("Elution 4 recapture", 60, 0, 0, NO_FLOW, buffer5_water, None),
        ("Elution 4 recapture loop", 200, 30, 0, LOOP, buffer5_water, None),
        ("Elution 4 pause", 154, 0, 0, NO_FLOW, buffer5_water, None),
        ("Elution 5 fill", 52, 40, 0, LINE, buffer3, frac_elution_4),
        ("Elution 5 resuspend", 90, 0, 60, NO_FLOW, buffer5_water, None),
        ("Elution 5 resuspend loop", 300, 30, 40, LOOP, buffer5_water, None),
        ("Elution 5 recapture", 60, 0, 0, NO_FLOW, buffer5_water, None),
        ("Elution 5 recapture loop", 200, 30, 0, LOOP, buffer5_water, None),
        ("Elution 5 pause", 170, 0, 0, NO_FLOW, buffer5_water, None),
        ("Regeneration", 52, 25, 0, LINE, buffer1, frac_elution_5),
    ]

    pump_percentages = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90], dtype=np.float64)
    upper_curve = np.array([0, 91, 220, 330, 465, 621, 783, 916, 1075, 1236, 1420, 1495, 1750, 1830, 1858, 1865, 1860, 1865, 1895], dtype=np.float64)
    lower_curve = np.array([0, 117, 260, 380, 505, 654, 783, 920, 1080, 1238, 1424, 1501, 1758, 1838, 1861, 1875, 1868, 1875, 1910], dtype=np.float64)
    flow_curve = (upper_curve + lower_curve) / 2.0 / 6e7

    mixing_percentages = np.array([0, 20, 40, 60, 80], dtype=np.float64)
    d_ax_curve = np.array([0.0, 1.3e-5, 5e-5, 6.9e-5, 8.4e-5], dtype=np.float64)

    cumulative_ends = []
    total_duration = 0.0
    for section in recipe:
        total_duration += section[1]
        cumulative_ends.append(total_duration)

    def section_idx(t: float):
        idx = bisect_right(cumulative_ends, t)
        return None if idx >= len(recipe) else idx

    def interp_pump(p: float) -> float:
        return float(np.interp(p, pump_percentages, flow_curve))

    def interp_dax(p: float) -> float:
        return float(np.interp(p, mixing_percentages, d_ax_curve))

    def flowRateFunction(t: float) -> float:
        idx = section_idx(t)
        if idx is None:
            return interp_pump(0.0)
        return interp_pump(recipe[idx][2])

    def pc_D_ax_function(t: float) -> float:
        idx = section_idx(t)
        if idx is None:
            return interp_dax(0.0)
        return interp_dax(recipe[idx][3])

    def capture_function(t: float) -> float:
        idx = section_idx(t)
        if idx is None:
            return 1.0
        section = recipe[idx]
        if section[3] == 0.0:
            return 1.0
        t_start = cumulative_ends[idx] - section[1]
        t_in_section = t - t_start
        uncapture_rate = math.log(2.0) / 5.0
        if t_in_section < 10.0:
            ramp = t_in_section / 10.0
            return -(ramp * ramp) * uncapture_rate
        return -uncapture_rate

    def inlet_function(t: float):
        idx = section_idx(t)
        return recipe[-1][5] if idx is None else recipe[idx][5]

    def flowRate_function(t: float) -> float:
        idx = section_idx(t)
        if idx is None:
            return interp_pump(0.0)
        cfg = recipe[idx][4]
        if cfg in (LINE, LOOP):
            return interp_pump(recipe[idx][2])
        return 0.0

    def flowRate_function_line(t: float) -> float:
        idx = section_idx(t)
        if idx is None:
            return interp_pump(0.0)
        return interp_pump(recipe[idx][2]) if recipe[idx][4] == LINE else 0.0

    def flowRate_function_loop(t: float) -> float:
        idx = section_idx(t)
        if idx is None:
            return interp_pump(0.0)
        return interp_pump(recipe[idx][2]) if recipe[idx][4] == LOOP else 0.0

    def flowRate_function_frac(t: float, fraction):
        idx = section_idx(t)
        if idx is None:
            return interp_pump(0.0)
        return interp_pump(recipe[idx][2]) if recipe[idx][6] is fraction else 0.0

    # =============================================================
    # ========================== UNIT OPERATIONS ==================
    # =============================================================
    inlet = fs3.Inlet(rs, inlet_function)

    d_pipes = 0.0096
    a_cross_pipes = math.pi / 4.0 * d_pipes * d_pipes
    l_inlet_pipe = 0.67 + (0.00005 / a_cross_pipes)
    scale_n_cells = lambda base: max(1, int(round(base * discretization_factor)))

    pipe_inlet = fs3.Pipe(rs, scale_n_cells(10), a_cross_pipes, l_inlet_pipe, flowRateFunction, 0.0)
    pipe_inlet.set_const_initial_concentration(buffer5_water)

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
    pc.set_const_initial_concentration(buffer5_water)

    v_dead = 4.113e-5 + 8.381e-5
    dead_volume = fs3.Volume(rs, v_dead)
    dead_volume.set_const_initial_concentration(buffer5_water)

    pipe_outlet = fs3.Pipe(rs, scale_n_cells(10), a_cross_pipes, 1.02, flowRateFunction, 0.0)
    pipe_outlet.set_const_initial_concentration(buffer5_water)

    pipe_loop = fs3.Pipe(rs, scale_n_cells(10), a_cross_pipes, 1.73, flowRate_function_loop, 0.0)
    pipe_loop.set_const_initial_concentration(buffer5_water)

    print(f"Total duration of the process: {total_duration} seconds.")

    process = fs3.Process(
        cs,
        [
            inlet,
            pipe_inlet,
            pc,
            dead_volume,
            pipe_outlet,
            pipe_loop,
            frac_feed_1,
            frac_feed_2,
            frac_wash_1,
            frac_wash_2,
            frac_wash_3,
            frac_elution_1,
            frac_elution_2,
            frac_elution_3,
            frac_elution_4,
            frac_elution_5,
        ],
    )

    process.add_connection(inlet.exit(), pipe_inlet.entry(), flowRate_function_line)
    process.add_connection(pipe_inlet.exit(), pc.entry(), flowRateFunction)
    process.add_connection(pc.exit(), dead_volume.entry(), flowRateFunction)
    process.add_connection(dead_volume.exit(), pipe_outlet.entry(), flowRateFunction)

    process.add_connection(pipe_outlet.exit(), frac_feed_1.entry(), lambda t: flowRate_function_frac(t, frac_feed_1))
    process.add_connection(pipe_outlet.exit(), frac_feed_2.entry(), lambda t: flowRate_function_frac(t, frac_feed_2))
    process.add_connection(pipe_outlet.exit(), frac_wash_1.entry(), lambda t: flowRate_function_frac(t, frac_wash_1))
    process.add_connection(pipe_outlet.exit(), frac_wash_2.entry(), lambda t: flowRate_function_frac(t, frac_wash_2))
    process.add_connection(pipe_outlet.exit(), frac_wash_3.entry(), lambda t: flowRate_function_frac(t, frac_wash_3))
    process.add_connection(pipe_outlet.exit(), frac_elution_1.entry(), lambda t: flowRate_function_frac(t, frac_elution_1))
    process.add_connection(pipe_outlet.exit(), frac_elution_2.entry(), lambda t: flowRate_function_frac(t, frac_elution_2))
    process.add_connection(pipe_outlet.exit(), frac_elution_3.entry(), lambda t: flowRate_function_frac(t, frac_elution_3))
    process.add_connection(pipe_outlet.exit(), frac_elution_4.entry(), lambda t: flowRate_function_frac(t, frac_elution_4))
    process.add_connection(pipe_outlet.exit(), frac_elution_5.entry(), lambda t: flowRate_function_frac(t, frac_elution_5))

    process.add_connection(pipe_outlet.exit(), pipe_loop.entry(), flowRate_function_loop)
    process.add_connection(pipe_loop.exit(), pipe_inlet.entry(), flowRate_function_loop)

    # =============================================================
    # ========================== OBSERVER =========================
    # =============================================================
    solver = fs3.Solver(process, solver_type)

    dt_obs = 2.0
    n_obs_time_steps = int(total_duration / dt_obs) + 1
    print(f"Number of observation time steps: {n_obs_time_steps}")

    compute_errors = True
    obs_pipe_inlet = fs3.TimeSeriesObserver(0.0, total_duration, n_obs_time_steps, pipe_inlet.all(), solver, save_obs, compute_errors)
    obs_pc_liquid = fs3.TimeSeriesObserver(0.0, total_duration, n_obs_time_steps, pc.liquid(), solver, save_obs, compute_errors)
    obs_pc_slurry = fs3.TimeSeriesObserver(0.0, total_duration, n_obs_time_steps, pc.slurry(), solver, save_obs, compute_errors)
    obs_pipe_outlet = fs3.TimeSeriesObserver(0.0, total_duration, n_obs_time_steps, pipe_outlet.all(), solver, save_obs, compute_errors)
    obs_pipe_loop = fs3.TimeSeriesObserver(0.0, total_duration, n_obs_time_steps, pipe_loop.all(), solver, save_obs, compute_errors)
    obs_frac_feed_1 = fs3.TimeSeriesObserver(0.0, total_duration, n_obs_time_steps, frac_feed_1.all(), solver, save_obs, compute_errors)
    obs_frac_feed_2 = fs3.TimeSeriesObserver(0.0, total_duration, n_obs_time_steps, frac_feed_2.all(), solver, save_obs, compute_errors)
    obs_frac_wash_1 = fs3.TimeSeriesObserver(0.0, total_duration, n_obs_time_steps, frac_wash_1.all(), solver, save_obs, compute_errors)
    obs_frac_wash_2 = fs3.TimeSeriesObserver(0.0, total_duration, n_obs_time_steps, frac_wash_2.all(), solver, save_obs, compute_errors)
    obs_frac_wash_3 = fs3.TimeSeriesObserver(0.0, total_duration, n_obs_time_steps, frac_wash_3.all(), solver, save_obs, compute_errors)
    obs_frac_elution_1 = fs3.TimeSeriesObserver(0.0, total_duration, n_obs_time_steps, frac_elution_1.all(), solver, save_obs, compute_errors)
    obs_frac_elution_2 = fs3.TimeSeriesObserver(0.0, total_duration, n_obs_time_steps, frac_elution_2.all(), solver, save_obs, compute_errors)
    obs_frac_elution_3 = fs3.TimeSeriesObserver(0.0, total_duration, n_obs_time_steps, frac_elution_3.all(), solver, save_obs, compute_errors)
    obs_frac_elution_4 = fs3.TimeSeriesObserver(0.0, total_duration, n_obs_time_steps, frac_elution_4.all(), solver, save_obs, compute_errors)
    obs_frac_elution_5 = fs3.TimeSeriesObserver(0.0, total_duration, n_obs_time_steps, frac_elution_5.all(), solver, save_obs, compute_errors)
    obs_pc_outlet = fs3.TimeSeriesObserver(0.0, total_duration, n_obs_time_steps, pc.exit(), solver, save_obs)

    # =============================================================
    # ========================== SOLVING ==========================
    # =============================================================
    solver.solve(total_duration, timeout_seconds)
    t_solve_duration = solver.get_solve_time()

    # =============================================================
    # ========================== RESULTS ==========================
    # =============================================================
    print("Max pH error: <not exposed in Python bindings>")
    internal_time_stamps = solver.get_internal_time_stamps()

    run_dir = Path.cwd() / datetime.now().strftime("run_%Y-%m-%d_%H-%M-%S")
    obs_dir = run_dir / "obs"
    if save_obs:
        (run_dir / "logs").mkdir(parents=True, exist_ok=True)
        (run_dir / "bench").mkdir(parents=True, exist_ok=True)
        obs_dir.mkdir(parents=True, exist_ok=True)

        np.save(obs_dir / "internal_timestamps.npy", np.asarray(internal_time_stamps, dtype=np.float64))

        obs_unit_ops = obs_dir / "unit_operations.npz"
        obs_pipe_inlet.save_to_npz(str(obs_unit_ops), "pipe_inlet", "a")
        obs_pc_liquid.save_to_npz(str(obs_unit_ops), "pc_liquid", "a")
        obs_pc_slurry.save_to_npz(str(obs_unit_ops), "pc_slurry", "a")
        obs_pipe_outlet.save_to_npz(str(obs_unit_ops), "pipe_outlet", "a")
        obs_pipe_loop.save_to_npz(str(obs_unit_ops), "pipe_loop", "a")
        obs_frac_feed_1.save_to_npz(str(obs_unit_ops), "frac_feed_1", "a")
        obs_frac_feed_2.save_to_npz(str(obs_unit_ops), "frac_feed_2", "a")
        obs_frac_wash_1.save_to_npz(str(obs_unit_ops), "frac_wash_1", "a")
        obs_frac_wash_2.save_to_npz(str(obs_unit_ops), "frac_wash_2", "a")
        obs_frac_wash_3.save_to_npz(str(obs_unit_ops), "frac_wash_3", "a")
        obs_frac_elution_1.save_to_npz(str(obs_unit_ops), "frac_elution_1", "a")
        obs_frac_elution_2.save_to_npz(str(obs_unit_ops), "frac_elution_2", "a")
        obs_frac_elution_3.save_to_npz(str(obs_unit_ops), "frac_elution_3", "a")
        obs_frac_elution_4.save_to_npz(str(obs_unit_ops), "frac_elution_4", "a")
        obs_frac_elution_5.save_to_npz(str(obs_unit_ops), "frac_elution_5", "a")

        obs_pc_outlet.save_to_npy(str(obs_dir / "pc_outlet.npy"))

        cell_index = pipe_outlet.n_cells // 2
        pH_activity = fs3.convert_to_pH(obs_pipe_outlet, cell_index, pipe_outlet.cell_volume, cs, activity_model)
        np.save(obs_dir / "pipe_outlet_middle_cell_pH_activity.npy", np.asarray(pH_activity, dtype=np.float64))

        pH_concentration = fs3.convert_to_pH(obs_pipe_outlet, cell_index, pipe_outlet.cell_volume, cs, fs3.NoActivityModel(cs))
        np.save(obs_dir / "pipe_outlet_middle_cell_pH_concentration.npy", np.asarray(pH_concentration, dtype=np.float64))

    print(f"Solved with {len(internal_time_stamps)} internal time stamps in {t_solve_duration} seconds.")
    if t_solve_duration > timeout_seconds:
        print(f"WARNING: Solver timed out after {timeout_seconds} seconds at t={solver.get_t()} / {total_duration}")

    return (float("nan"), t_solve_duration, len(internal_time_stamps))


if __name__ == "__main__":
    run_HGMS_hIgG_separation(
        1e3,
        0.1,
        True,
        float("inf"),
        fs3.SolverType.ERK,
        0.1,
    )

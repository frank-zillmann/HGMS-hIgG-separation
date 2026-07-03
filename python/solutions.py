"""Inlet solutions: build and equilibrate the feed and buffers."""

import math
from typing import Dict

import numpy as np

import fs3
from components import mnp_ns, mnp_specific_surface
from reactions import langmuir_q_max_from_pH


def _equilibrate(cs, rs, concentrations, t_duration, solver_type, timeout_seconds, label, expected):
    concentrations, error = fs3.one_cell_reaction(rs, concentrations, t_duration, solver_type, timeout_seconds)
    print(f"pH {label}: {-math.log10(concentrations[cs.get_idx('H⁺')] * 1e-3)} (should be ~{expected})")
    print(f"{label} reaction error: {error}")
    return concentrations


def build_solutions(cs, rs, solver_type, timeout_seconds=float("inf")) -> Dict[str, np.ndarray]:
    """Build and equilibrate the inlet solutions; returns ``{name: concentrations}``."""
    c_w = 1000.0 / cs["H₂O"].molar_mass
    t_reaction_duration = 100.0

    def blank():
        return np.zeros(cs.n_components, dtype=np.float64)

    # ----- Feed -----
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

    feed = blank()
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

    feed = _equilibrate(cs, rs, feed, t_reaction_duration * 10.0, solver_type, timeout_seconds, "feed", "7.3")
    pH_feed = -math.log10(feed[cs.get_idx("H⁺")] * 1e-3)
    q_feed = feed[cs.get_idx("MNP-hIgG")] / feed[cs.get_idx("MNP1") : cs.get_idx("MNP10") + 1].sum()
    print(f"q feed: {q_feed} q_max at this pH: {langmuir_q_max_from_pH(pH_feed)}")

    # ----- Buffer 1 -----
    buffer1 = blank()
    buffer1[cs.get_idx("H₂O")] = c_w
    buffer1[cs.get_idx("H⁺")] = 0.019156e3
    buffer1[cs.get_idx("Na⁺")] = 0.15e3
    buffer1[cs.get_idx("Cl⁻")] = 0.019156e3 + 0.15e3
    buffer1[cs.get_idx("Tris")] = 0.02e3
    buffer1 = _equilibrate(cs, rs, buffer1, t_reaction_duration, solver_type, timeout_seconds, "B1", "7")

    # ----- Buffer 2 -----
    buffer2 = blank()
    buffer2[cs.get_idx("H₂O")] = c_w
    buffer2[cs.get_idx("H⁺")] = 0.072107e3
    buffer2[cs.get_idx("Na⁺")] = 0.2e3
    buffer2[cs.get_idx("Cl⁻")] = 0.072107e3
    buffer2[cs.get_idx("Ac⁻")] = 0.2e3
    buffer2 = _equilibrate(cs, rs, buffer2, t_reaction_duration, solver_type, timeout_seconds, "B2", "4.75")

    # ----- Buffer 3 -----
    buffer3 = blank()
    buffer3[cs.get_idx("H₂O")] = c_w
    buffer3[cs.get_idx("H⁺")] = 0.20221e3
    buffer3[cs.get_idx("Na⁺")] = 0.2e3
    buffer3[cs.get_idx("Cl⁻")] = 0.20221e3
    buffer3[cs.get_idx("Ac⁻")] = 0.2e3
    buffer3 = _equilibrate(cs, rs, buffer3, t_reaction_duration, solver_type, timeout_seconds, "B3", "2.5")

    # ----- Water (B5) -----
    water = blank()
    water[cs.get_idx("H₂O")] = c_w
    water = _equilibrate(cs, rs, water, t_reaction_duration, solver_type, timeout_seconds, "B5", "7.0 / 6.5")

    return {
        "Feed": feed,
        "Buffer 1": buffer1,
        "Buffer 2": buffer2,
        "Buffer 3": buffer3,
        "Water (B5)": water,
    }

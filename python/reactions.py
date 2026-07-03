"""Reaction system (acid/base equilibria, Gouy-Chapman surface, Langmuir binding)."""

import numpy as np

import fs3
from components import mnp_ns

# Langmuir isotherm support points: pH -> binding constant / capacity.
langmuir_ph = np.array([2.0, 2.5, 2.6, 3.5, 3.7, 4.5, 7.0], dtype=np.float64)
langmuir_k_b = np.array([0.03, 2.14, 2.28, 3.22, 3.84, 35.15, 36.52], dtype=np.float64)
langmuir_q_max = np.array([0.00, 0.01, 0.02, 0.05, 0.13, 0.28, 0.31], dtype=np.float64)


def langmuir_k_b_from_pH(pH: float) -> float:
    return float(np.interp(pH, langmuir_ph, langmuir_k_b))


def langmuir_q_max_from_pH(pH: float) -> float:
    return float(np.interp(pH, langmuir_ph, langmuir_q_max))


def build_reaction_system(cs: fs3.ComponentSystem, tau_reaction: float):
    """Build the reaction system; returns ``(reaction_system, activity_model)``."""
    activity_model = fs3.TruesdellJonesActivityModel(cs)
    rs = fs3.ReactionSystem(cs, activity_model)

    c_w = 1000.0 / cs["H₂O"].molar_mass
    ks_w = (10.0 ** (-14) * 1e6) / c_w
    ks_tris_h_plus = 10.0 ** (-8.3) * 1e3
    ks_ac_h = 10.0 ** (-4.75) * 1e3
    ks_gly_h2_plus = 10.0 ** (-2.35) * 1e3
    ks_gly_h = 10.0 ** (-9.78) * 1e3
    ks_mnp_oh2_plus = 10.0 ** (-(7.9 - 0.5 * 2.5)) * 1e3
    ks_mnp_oh = 10.0 ** (-(7.9 + 0.5 * 2.5)) * 1e3

    h_plus = cs.get_idx("H⁺")
    rs.add(fs3.mass_action_law_inverse_rate_prediction(cs, "H₂O <=> H⁺ + OH⁻", ks_w, tau_reaction, h_plus))
    rs.add(fs3.mass_action_law_inverse_rate_prediction(cs, "TrisH⁺ <=> Tris + H⁺", ks_tris_h_plus, tau_reaction, h_plus))
    rs.add(fs3.mass_action_law_inverse_rate_prediction(cs, "AcH <=> Ac⁻ + H⁺", ks_ac_h, tau_reaction, h_plus))
    rs.add(fs3.mass_action_law_inverse_rate_prediction(cs, "GlyH₂⁺ <=> GlyH + H⁺", ks_gly_h2_plus, tau_reaction, h_plus))
    rs.add(fs3.mass_action_law_inverse_rate_prediction(cs, "GlyH <=> Gly⁻ + H⁺", ks_gly_h, tau_reaction, h_plus))

    gouy_chapman_model = fs3.GouyChapmanModel(cs, "MNP-OH₂⁺", "MNP-OH", "MNP-O⁻", "H⁺", mnp_ns)
    reaction_oh2_plus = fs3.mass_action_law_inverse_rate_prediction(cs, "MNP-OH₂⁺ <=> MNP-OH + H⁺", ks_mnp_oh2_plus, tau_reaction, h_plus)
    reaction_oh = fs3.mass_action_law_inverse_rate_prediction(cs, "MNP-OH <=> MNP-O⁻ + H⁺ ", ks_mnp_oh, tau_reaction, h_plus)
    rs.add(fs3.wrap_reaction(reaction_oh2_plus, gouy_chapman_model))
    rs.add(fs3.wrap_reaction(reaction_oh, gouy_chapman_model))

    rs.add(
        fs3.langmuir_binding_reaction(
            h_plus,
            cs.get_idx("MNP1"),
            cs.get_idx("MNP10"),
            cs.get_idx("MNP-hIgG"),
            cs.get_idx("hIgG"),
            cs["hIgG"].molar_mass,
            langmuir_k_b_from_pH,
            langmuir_q_max_from_pH,
            tau_reaction,
        )
    )
    return rs, activity_model

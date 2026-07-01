"""Shared FS³ model: component system, reactions, solutions and the recipe schema.

Everything that the standalone simulation script (``HGMS_hIgG_separation.py``), the UI and the
plot scripts need to agree on lives here, so it is defined in exactly one place.
"""

import math
from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Optional

import numpy as np

import fs3

# Surface-site density and specific surface area of the MNPs.
mnp_ns = 1.08e-5
mnp_specific_surface = 1e5

# Langmuir isotherm support points: pH -> binding constant / capacity.
langmuir_ph = np.array([2.0, 2.5, 2.6, 3.5, 3.7, 4.5, 7.0], dtype=np.float64)
langmuir_k_b = np.array([0.03, 2.14, 2.28, 3.22, 3.84, 35.15, 36.52], dtype=np.float64)
langmuir_q_max = np.array([0.00, 0.01, 0.02, 0.05, 0.13, 0.28, 0.31], dtype=np.float64)


def langmuir_k_b_from_pH(pH: float) -> float:
    return float(np.interp(pH, langmuir_ph, langmuir_k_b))


def langmuir_q_max_from_pH(pH: float) -> float:
    return float(np.interp(pH, langmuir_ph, langmuir_q_max))


def build_component_system() -> fs3.ComponentSystem:
    return fs3.ComponentSystem(
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


# =============================================================
# ========================== RECIPE SCHEMA ====================
# =============================================================
class FlowSheetConfig(Enum):
    NO_FLOW = 0
    LINE = 1
    LOOP = 2


@dataclass(slots=True)
class RecipeStep:
    name: str
    t_duration: float
    pump_percentage: float
    mixing_percentage: float
    flow_sheet_configuration: FlowSheetConfig
    inlet: str  # solution name, resolved against build_solutions() output
    fraction: Optional[str]


FRACTIONS = [
    "Feed 1", "Feed 2",
    "Wash 1", "Wash 2", "Wash 3",
    "Elution 1", "Elution 2", "Elution 3", "Elution 4", "Elution 5",
]

FLOW_SHEET_LABELS: Dict[FlowSheetConfig, str] = {
    FlowSheetConfig.NO_FLOW: "No flow",
    FlowSheetConfig.LINE: "Line",
    FlowSheetConfig.LOOP: "Loop",
}


def default_recipe() -> List[RecipeStep]:
    return [
        RecipeStep("Load 1", 1199, 20, 0, FlowSheetConfig.LINE, "Feed", "Feed 1"),
        RecipeStep("Load 1 pause", 34, 0, 0, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        RecipeStep("Load 2", 96, 20, 0, FlowSheetConfig.LINE, "Feed", "Feed 1"),
        RecipeStep("Load 2 pause", 33, 0, 0, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        RecipeStep("Wash 1 fill", 50, 40, 0, FlowSheetConfig.LINE, "Buffer 1", "Feed 2"),
        RecipeStep("Wash 1 resuspend", 65, 0, 60, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        RecipeStep("Wash 1 resuspend loop", 51, 40, 50, FlowSheetConfig.LOOP, "Water (B5)", None),
        RecipeStep("Wash 1 recapture", 69, 0, 0, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        RecipeStep("Wash 1 recapture loop", 196, 30, 0, FlowSheetConfig.LOOP, "Water (B5)", None),
        RecipeStep("Wash 1 pause", 6, 0, 0, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        RecipeStep("Wash 2 fill", 50, 40, 0, FlowSheetConfig.LINE, "Buffer 1", "Wash 1"),
        RecipeStep("Wash 2 resuspend", 65, 0, 60, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        RecipeStep("Wash 2 resuspend loop", 51, 40, 50, FlowSheetConfig.LOOP, "Water (B5)", None),
        RecipeStep("Wash 2 recapture", 69, 0, 0, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        RecipeStep("Wash 2 recapture loop", 196, 30, 0, FlowSheetConfig.LOOP, "Water (B5)", None),
        RecipeStep("Wash 2 pause", 14, 0, 0, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        RecipeStep("Wash 3 fill", 50, 40, 0, FlowSheetConfig.LINE, "Buffer 2", "Wash 2"),
        RecipeStep("Wash 3 resuspend", 65, 0, 60, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        RecipeStep("Wash 3 resuspend loop", 51, 40, 50, FlowSheetConfig.LOOP, "Water (B5)", None),
        RecipeStep("Wash 3 recapture", 69, 0, 0, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        RecipeStep("Wash 3 recapture loop", 196, 30, 0, FlowSheetConfig.LOOP, "Water (B5)", None),
        RecipeStep("Wash 3 pause", 13, 0, 0, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        RecipeStep("Elution 1 fill", 52, 40, 0, FlowSheetConfig.LINE, "Buffer 3", "Wash 3"),
        RecipeStep("Elution 1 resuspend", 90, 0, 60, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        RecipeStep("Elution 1 resuspend loop", 300, 30, 40, FlowSheetConfig.LOOP, "Water (B5)", None),
        RecipeStep("Elution 1 recapture", 60, 0, 0, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        RecipeStep("Elution 1 recapture loop", 200, 30, 0, FlowSheetConfig.LOOP, "Water (B5)", None),
        RecipeStep("Elution 1 pause", 7, 0, 0, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        RecipeStep("Elution 2 fill", 52, 40, 0, FlowSheetConfig.LINE, "Buffer 3", "Elution 1"),
        RecipeStep("Elution 2 resuspend", 90, 0, 60, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        RecipeStep("Elution 2 resuspend loop", 300, 30, 40, FlowSheetConfig.LOOP, "Water (B5)", None),
        RecipeStep("Elution 2 recapture", 60, 0, 0, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        RecipeStep("Elution 2 recapture loop", 200, 30, 0, FlowSheetConfig.LOOP, "Water (B5)", None),
        RecipeStep("Elution 2 pause", 9, 0, 0, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        RecipeStep("Elution 3 fill", 52, 40, 0, FlowSheetConfig.LINE, "Buffer 3", "Elution 2"),
        RecipeStep("Elution 3 resuspend", 90, 0, 60, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        RecipeStep("Elution 3 resuspend loop", 300, 30, 40, FlowSheetConfig.LOOP, "Water (B5)", None),
        RecipeStep("Elution 3 recapture", 60, 0, 0, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        RecipeStep("Elution 3 recapture loop", 200, 30, 0, FlowSheetConfig.LOOP, "Water (B5)", None),
        RecipeStep("Elution 3 pause", 10, 0, 0, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        RecipeStep("Elution 4 fill", 52, 40, 0, FlowSheetConfig.LINE, "Buffer 3", "Elution 3"),
        RecipeStep("Elution 4 resuspend", 90, 0, 60, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        RecipeStep("Elution 4 resuspend loop", 300, 30, 40, FlowSheetConfig.LOOP, "Water (B5)", None),
        RecipeStep("Elution 4 recapture", 60, 0, 0, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        RecipeStep("Elution 4 recapture loop", 200, 30, 0, FlowSheetConfig.LOOP, "Water (B5)", None),
        RecipeStep("Elution 4 pause", 154, 0, 0, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        RecipeStep("Elution 5 fill", 52, 40, 0, FlowSheetConfig.LINE, "Buffer 3", "Elution 4"),
        RecipeStep("Elution 5 resuspend", 90, 0, 60, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        RecipeStep("Elution 5 resuspend loop", 300, 30, 40, FlowSheetConfig.LOOP, "Water (B5)", None),
        RecipeStep("Elution 5 recapture", 60, 0, 0, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        RecipeStep("Elution 5 recapture loop", 200, 30, 0, FlowSheetConfig.LOOP, "Water (B5)", None),
        RecipeStep("Elution 5 pause", 170, 0, 0, FlowSheetConfig.NO_FLOW, "Water (B5)", None),
        RecipeStep("Regeneration", 52, 25, 0, FlowSheetConfig.LINE, "Buffer 1", "Elution 5"),
    ]


def total_duration(recipe: List[RecipeStep]) -> float:
    return float(sum(step.t_duration for step in recipe))


def observation_time_stamps(recipe: List[RecipeStep], dt_obs: float = 2.0) -> np.ndarray:
    return np.arange(0.0, total_duration(recipe) + 1e-5, dt_obs)


# Default run directory used by the plot/animation scripts when run standalone:
# <DEFAULT_RUN>/obs holds the simulation output, <DEFAULT_RUN>/plots the figures.
DEFAULT_RUN = "data/run_fitted_no_MNP_hydroxyl_reactions_ERK_for_time_steps"
path_to_obs = f"{DEFAULT_RUN}/obs"
path_to_save = f"{DEFAULT_RUN}/plots"

# Runs to overlay in plot_fractions.py, as {label: obs directory}.
path_to_obs_dict = {
    "FS³, fitted Langmuir, no MNP hydroxyl reactions": "data/run_fitted_no_MNP_hydroxyl_reactions/obs",
    "FS³, fitted Langmuir, no MNP hydroxyl reactions, tau=1, discretization factor = 0.5": "data/run_fitted_no_MNP_hydroxyl_reactions_tau_1_discretization_factor_0,5_ERK/obs",
    "FS³, fitted Langmuir, no MNP hydroxyl reactions, tau=5, discretization factor = 0.2": "data/run_fitted_no_MNP_hydroxyl_reactions_tau_5_discretization_factor_0,2_ERK/obs",
}

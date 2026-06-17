from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Optional

import numpy as np


class FlowSheetConfig(Enum):
    NO_FLOW = 0
    LINE = 1
    LOOP = 2


@dataclass(slots=True)
class Solution:
    """A named inlet solution and its component concentrations."""
    name: str
    concentrations: np.ndarray


@dataclass(slots=True)
class RecipeStep:
    name: str
    t_duration: float
    pump_percentage: float
    mixing_percentage: float
    flow_sheet_configuration: FlowSheetConfig
    inlet: Solution
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


def default_recipe(solutions: Dict[str, Solution]) -> List[RecipeStep]:
    return [
        RecipeStep("Load 1", 1199, 20, 0, FlowSheetConfig.LINE, solutions["Feed"], "Feed 1"),
        RecipeStep("Load 1 pause", 34, 0, 0, FlowSheetConfig.NO_FLOW, solutions["Water (B5)"], None),
        RecipeStep("Load 2", 96, 20, 0, FlowSheetConfig.LINE, solutions["Feed"], "Feed 1"),
        RecipeStep("Load 2 pause", 33, 0, 0, FlowSheetConfig.NO_FLOW, solutions["Water (B5)"], None),
        RecipeStep("Wash 1 fill", 50, 40, 0, FlowSheetConfig.LINE, solutions["Buffer 1"], "Feed 2"),
        RecipeStep("Wash 1 resuspend", 65, 0, 60, FlowSheetConfig.NO_FLOW, solutions["Water (B5)"], None),
        RecipeStep("Wash 1 resuspend loop", 51, 40, 50, FlowSheetConfig.LOOP, solutions["Water (B5)"], None),
        RecipeStep("Wash 1 recapture", 69, 0, 0, FlowSheetConfig.NO_FLOW, solutions["Water (B5)"], None),
        RecipeStep("Wash 1 recapture loop", 196, 30, 0, FlowSheetConfig.LOOP, solutions["Water (B5)"], None),
        RecipeStep("Wash 1 pause", 6, 0, 0, FlowSheetConfig.NO_FLOW, solutions["Water (B5)"], None),
        RecipeStep("Wash 2 fill", 50, 40, 0, FlowSheetConfig.LINE, solutions["Buffer 1"], "Wash 1"),
        RecipeStep("Wash 2 resuspend", 65, 0, 60, FlowSheetConfig.NO_FLOW, solutions["Water (B5)"], None),
        RecipeStep("Wash 2 resuspend loop", 51, 40, 50, FlowSheetConfig.LOOP, solutions["Water (B5)"], None),
        RecipeStep("Wash 2 recapture", 69, 0, 0, FlowSheetConfig.NO_FLOW, solutions["Water (B5)"], None),
        RecipeStep("Wash 2 recapture loop", 196, 30, 0, FlowSheetConfig.LOOP, solutions["Water (B5)"], None),
        RecipeStep("Wash 2 pause", 14, 0, 0, FlowSheetConfig.NO_FLOW, solutions["Water (B5)"], None),
        RecipeStep("Wash 3 fill", 50, 40, 0, FlowSheetConfig.LINE, solutions["Buffer 2"], "Wash 2"),
        RecipeStep("Wash 3 resuspend", 65, 0, 60, FlowSheetConfig.NO_FLOW, solutions["Water (B5)"], None),
        RecipeStep("Wash 3 resuspend loop", 51, 40, 50, FlowSheetConfig.LOOP, solutions["Water (B5)"], None),
        RecipeStep("Wash 3 recapture", 69, 0, 0, FlowSheetConfig.NO_FLOW, solutions["Water (B5)"], None),
        RecipeStep("Wash 3 recapture loop", 196, 30, 0, FlowSheetConfig.LOOP, solutions["Water (B5)"], None),
        RecipeStep("Wash 3 pause", 13, 0, 0, FlowSheetConfig.NO_FLOW, solutions["Water (B5)"], None),
        RecipeStep("Elution 1 fill", 52, 40, 0, FlowSheetConfig.LINE, solutions["Buffer 3"], "Wash 3"),
        RecipeStep("Elution 1 resuspend", 90, 0, 60, FlowSheetConfig.NO_FLOW, solutions["Water (B5)"], None),
        RecipeStep("Elution 1 resuspend loop", 300, 30, 40, FlowSheetConfig.LOOP, solutions["Water (B5)"], None),
        RecipeStep("Elution 1 recapture", 60, 0, 0, FlowSheetConfig.NO_FLOW, solutions["Water (B5)"], None),
        RecipeStep("Elution 1 recapture loop", 200, 30, 0, FlowSheetConfig.LOOP, solutions["Water (B5)"], None),
        RecipeStep("Elution 1 pause", 7, 0, 0, FlowSheetConfig.NO_FLOW, solutions["Water (B5)"], None),
        RecipeStep("Elution 2 fill", 52, 40, 0, FlowSheetConfig.LINE, solutions["Buffer 3"], "Elution 1"),
        RecipeStep("Elution 2 resuspend", 90, 0, 60, FlowSheetConfig.NO_FLOW, solutions["Water (B5)"], None),
        RecipeStep("Elution 2 resuspend loop", 300, 30, 40, FlowSheetConfig.LOOP, solutions["Water (B5)"], None),
        RecipeStep("Elution 2 recapture", 60, 0, 0, FlowSheetConfig.NO_FLOW, solutions["Water (B5)"], None),
        RecipeStep("Elution 2 recapture loop", 200, 30, 0, FlowSheetConfig.LOOP, solutions["Water (B5)"], None),
        RecipeStep("Elution 2 pause", 9, 0, 0, FlowSheetConfig.NO_FLOW, solutions["Water (B5)"], None),
        RecipeStep("Elution 3 fill", 52, 40, 0, FlowSheetConfig.LINE, solutions["Buffer 3"], "Elution 2"),
        RecipeStep("Elution 3 resuspend", 90, 0, 60, FlowSheetConfig.NO_FLOW, solutions["Water (B5)"], None),
        RecipeStep("Elution 3 resuspend loop", 300, 30, 40, FlowSheetConfig.LOOP, solutions["Water (B5)"], None),
        RecipeStep("Elution 3 recapture", 60, 0, 0, FlowSheetConfig.NO_FLOW, solutions["Water (B5)"], None),
        RecipeStep("Elution 3 recapture loop", 200, 30, 0, FlowSheetConfig.LOOP, solutions["Water (B5)"], None),
        RecipeStep("Elution 3 pause", 10, 0, 0, FlowSheetConfig.NO_FLOW, solutions["Water (B5)"], None),
        RecipeStep("Elution 4 fill", 52, 40, 0, FlowSheetConfig.LINE, solutions["Buffer 3"], "Elution 3"),
        RecipeStep("Elution 4 resuspend", 90, 0, 60, FlowSheetConfig.NO_FLOW, solutions["Water (B5)"], None),
        RecipeStep("Elution 4 resuspend loop", 300, 30, 40, FlowSheetConfig.LOOP, solutions["Water (B5)"], None),
        RecipeStep("Elution 4 recapture", 60, 0, 0, FlowSheetConfig.NO_FLOW, solutions["Water (B5)"], None),
        RecipeStep("Elution 4 recapture loop", 200, 30, 0, FlowSheetConfig.LOOP, solutions["Water (B5)"], None),
        RecipeStep("Elution 4 pause", 154, 0, 0, FlowSheetConfig.NO_FLOW, solutions["Water (B5)"], None),
        RecipeStep("Elution 5 fill", 52, 40, 0, FlowSheetConfig.LINE, solutions["Buffer 3"], "Elution 4"),
        RecipeStep("Elution 5 resuspend", 90, 0, 60, FlowSheetConfig.NO_FLOW, solutions["Water (B5)"], None),
        RecipeStep("Elution 5 resuspend loop", 300, 30, 40, FlowSheetConfig.LOOP, solutions["Water (B5)"], None),
        RecipeStep("Elution 5 recapture", 60, 0, 0, FlowSheetConfig.NO_FLOW, solutions["Water (B5)"], None),
        RecipeStep("Elution 5 recapture loop", 200, 30, 0, FlowSheetConfig.LOOP, solutions["Water (B5)"], None),
        RecipeStep("Elution 5 pause", 170, 0, 0, FlowSheetConfig.NO_FLOW, solutions["Water (B5)"], None),
        RecipeStep("Regeneration", 52, 25, 0, FlowSheetConfig.LINE, solutions["Buffer 1"], "Elution 5"),
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
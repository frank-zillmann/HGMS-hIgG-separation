"""Recipe schema: the editable process steps that drive a simulation."""

from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Optional


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


def total_duration(recipe: List[RecipeStep]) -> float:
    return float(sum(step.t_duration for step in recipe))


def phase_transitions(recipe: List[RecipeStep]) -> List[float]:
    """Cumulative end time of every step -> vertical guide lines matching the recipe."""
    ends, total = [], 0.0
    for step in recipe:
        total += step.t_duration
        ends.append(total)
    return ends


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

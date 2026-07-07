"""Defines the HGMS hIgG separation experiment.

This module only *defines* how to build and run an experiment; it never runs one on import
(no side effects, no threading). Callers -- ``run.py`` for batch use, ``ui.py`` for the live
UI -- drive an ``Experiment`` and read observations straight off it as numpy arrays.
"""

from dataclasses import dataclass, field
from typing import Dict, Optional, Sequence

import numpy as np

import fs3
from components import build_component_system
from observers import build_observers, fraction_masses, outlet_pH
from reactions import build_reaction_system
from recipe import RecipeStep, default_recipe, total_duration
from solutions import build_solutions
from unit_operations import build_process


@dataclass
class Experiment:
    solver: fs3.Solver
    observers: Dict[str, fs3.TimeSeriesObserver]
    unit_operations: Dict[str, object]
    cs: fs3.ComponentSystem
    activity_model: fs3.ActivityModelBase
    recipe: Sequence[RecipeStep]
    total_duration: float
    n_obs_time_steps: int
    timeout_seconds: float = float("inf")
    _snapshot_times: np.ndarray = field(init=False)

    def __post_init__(self):
        self._snapshot_times = np.linspace(0.0, self.total_duration, self.n_obs_time_steps)

    # ----- running -----
    def solve(self) -> None:
        self.solver.solve(self.total_duration, self.timeout_seconds)

    @property
    def t(self) -> float:
        return self.solver.get_t()

    @property
    def progress(self) -> float:
        return min(self.t / self.total_duration, 1.0)

    def n_filled(self) -> int:
        """Number of observation snapshots the solver has already reached."""
        return max(0, int(np.searchsorted(self._snapshot_times, self.t, side="right")) - 1)

    # ----- reading observations (one method per result) -----
    def times(self) -> np.ndarray:
        return self._snapshot_times

    def outlet_pH(self, activity: bool = True) -> np.ndarray:
        model = self.activity_model if activity else fs3.NoActivityModel(self.cs)
        return outlet_pH(self.observers, self.unit_operations["pipe_outlet"], self.cs, model)

    def live_pH(self):
        """Activity- and concentration-based pH of the reached snapshots (for live plotting)."""
        n = self.n_filled()
        times = self._snapshot_times[:n]
        ph_activity = self.outlet_pH(activity=True)[:n]
        ph_concentration = self.outlet_pH(activity=False)[:n]
        finite = np.isfinite(ph_activity)
        return times[finite], ph_activity[finite], ph_concentration[finite]

    def fraction_masses(self, component: str = "hIgG") -> Dict[str, float]:
        return fraction_masses(self.observers, self.cs.get_idx(component))

    def step_sizes(self):
        ts = np.asarray(self.solver.get_internal_time_stamps(), dtype=np.float64)
        return ts[:-1], np.diff(ts)


def build_experiment(
    recipe: Optional[Sequence[RecipeStep]] = None,
    *,
    discretization_factor: float = 0.1,
    solver_type=fs3.SolverType.ERK,
    tau_reaction: float = 0.1,
    dt_obs: float = 2.0,
    timeout_seconds: float = float("inf"),
    solutions: Optional[Dict[str, np.ndarray]] = None,
) -> Experiment:
    """Assemble a ready-to-solve experiment. Pass pre-equilibrated ``solutions`` (plain
    concentration arrays) to skip the ~2 s equilibration; the rest is always built fresh."""
    recipe = list(recipe) if recipe is not None else default_recipe()
    cs = build_component_system()
    rs, activity_model = build_reaction_system(cs, tau_reaction)
    if solutions is None:
        solutions = build_solutions(cs, rs, solver_type, timeout_seconds)

    process, unit_operations = build_process(rs, cs, recipe, solutions, discretization_factor)
    solver = fs3.Solver(process, solver_type)

    duration = total_duration(recipe)
    n_obs_time_steps = int(duration / dt_obs) + 1
    observers = build_observers(solver, unit_operations, duration, n_obs_time_steps)

    return Experiment(solver, observers, unit_operations, cs, activity_model, recipe, duration, n_obs_time_steps, timeout_seconds)


def run_HGMS_hIgG_separation(recipe: Optional[Sequence[RecipeStep]] = None, **kwargs) -> Experiment:
    """Build an experiment and solve it (blocking); returns the finished ``Experiment``."""
    experiment = build_experiment(recipe, **kwargs)
    experiment.solve()
    return experiment

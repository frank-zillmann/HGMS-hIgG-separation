"""Defines the HGMS hIgG separation experiment.

This module only *defines* how to build and run a simulation; it never runs one on import
(no side effects, no threading). Callers -- ``run.py`` for batch use, ``ui.py`` for the live
UI -- drive ``Simulation`` and read observations straight out of it as numpy arrays.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np

import fs3
from components import build_component_system
from observers import build_observers, fraction_masses, outlet_pH
from reactions import build_reaction_system
from recipe import RecipeStep, default_recipe, phase_transitions, total_duration
from solutions import build_solutions
from unit_operations import build_plant


@dataclass
class Results:
    """Plot-ready observations extracted from a finished simulation (no FS³ objects)."""
    times: np.ndarray
    ph_activity: np.ndarray
    ph_concentration: np.ndarray
    fraction_masses: Dict[str, float]
    step_times: np.ndarray
    step_sizes: np.ndarray
    transitions: List[float]


@dataclass
class Simulation:
    solver: fs3.Solver
    observers: Dict[str, fs3.TimeSeriesObserver]
    cs: fs3.ComponentSystem
    activity_model: fs3.ActivityModelBase
    plant: object
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

    # ----- reading observations -----
    def outlet_pH(self, activity: bool = True) -> np.ndarray:
        model = self.activity_model if activity else fs3.NoActivityModel(self.cs)
        return outlet_pH(self.observers, self.plant.pipe_outlet, self.cs, model)

    def live_pH(self) -> tuple:
        """Activity-based pH of the reached snapshots only (for live plotting during a run)."""
        n = self.n_filled()
        times, ph = self._snapshot_times[:n], self.outlet_pH(activity=True)[:n]
        finite = np.isfinite(ph)
        return times[finite], ph[finite]

    def fraction_masses(self, component: str = "hIgG") -> Dict[str, float]:
        return fraction_masses(self.observers, self.cs.get_idx(component))

    def step_sizes(self) -> tuple:
        ts = np.asarray(self.solver.get_internal_time_stamps(), dtype=np.float64)
        return ts[:-1], np.diff(ts)

    def results(self, component: str = "hIgG") -> Results:
        step_times, step_sizes = self.step_sizes()
        return Results(
            times=self._snapshot_times,
            ph_activity=self.outlet_pH(activity=True),
            ph_concentration=self.outlet_pH(activity=False),
            fraction_masses=self.fraction_masses(component),
            step_times=step_times,
            step_sizes=step_sizes,
            transitions=phase_transitions(self.recipe),
        )


def build_simulation(
    recipe: Optional[Sequence[RecipeStep]] = None,
    *,
    discretization_factor: float = 0.1,
    solver_type=fs3.SolverType.ERK,
    tau_reaction: float = 0.1,
    dt_obs: float = 2.0,
    timeout_seconds: float = float("inf"),
    solutions: Optional[Dict[str, np.ndarray]] = None,
) -> Simulation:
    """Assemble a ready-to-solve simulation. Pass pre-equilibrated ``solutions`` (plain
    concentration arrays) to skip the ~2 s equilibration; the rest is always built fresh."""
    recipe = list(recipe) if recipe is not None else default_recipe()
    cs = build_component_system()
    rs, activity_model = build_reaction_system(cs, tau_reaction)
    if solutions is None:
        solutions = build_solutions(cs, rs, solver_type, timeout_seconds)

    plant = build_plant(rs, cs, recipe, solutions, discretization_factor)
    solver = fs3.Solver(plant.process, solver_type)

    duration = total_duration(recipe)
    n_obs_time_steps = int(duration / dt_obs) + 1
    observers = build_observers(solver, plant, duration, n_obs_time_steps)

    return Simulation(solver, observers, cs, activity_model, plant, recipe, duration, n_obs_time_steps, timeout_seconds)


def run_HGMS_hIgG_separation(recipe: Optional[Sequence[RecipeStep]] = None, **kwargs) -> Simulation:
    """Build a simulation and solve it (blocking); returns the finished ``Simulation``."""
    sim = build_simulation(recipe, **kwargs)
    sim.solve()
    return sim

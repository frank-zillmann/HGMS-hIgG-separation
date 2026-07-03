"""Time-series observers and helpers to read observation data back out as numpy arrays."""

from pathlib import Path
from typing import Dict, Union

import numpy as np

import fs3
from recipe import FRACTIONS
from unit_operations import Plant


def build_observers(solver, plant: Plant, total_duration: float, n_obs_time_steps: int) -> Dict[str, fs3.TimeSeriesObserver]:
    """One observer per observed unit operation, registered with the solver (records in place)."""
    targets = {
        "pipe_inlet": plant.pipe_inlet.all(),
        "pc_liquid": plant.pc.liquid(),
        "pc_slurry": plant.pc.slurry(),
        "pipe_outlet": plant.pipe_outlet.all(),
        "pipe_loop": plant.pipe_loop.all(),
        **{name: frac.all() for name, frac in plant.fractions.items()},
    }
    return {
        name: fs3.TimeSeriesObserver(0.0, total_duration, n_obs_time_steps, target, solver, True, True)
        for name, target in targets.items()
    }


def observer_array(observer: fs3.TimeSeriesObserver) -> np.ndarray:
    """Full recorded series as ``(n_snapshots, n_cells, n_components)``."""
    values, shape = observer.get_all_snapshots_as_vector()
    return np.asarray(values, dtype=np.float64).reshape(shape)


def outlet_pH(observers, pipe_outlet, cs, activity_model) -> np.ndarray:
    cell_index = pipe_outlet.n_cells // 2
    return np.asarray(fs3.convert_to_pH(observers["pipe_outlet"], cell_index, pipe_outlet.cell_volume, cs, activity_model))


def fraction_masses(observers, component_idx: int) -> Dict[str, float]:
    """Captured mass [g] of one component in each fraction at the final time point."""
    return {name: float(observer_array(observers[name])[-1, 0, component_idx]) * 1000.0 for name in FRACTIONS}


def save_observations(sim, obs_dir: Union[str, Path]) -> Path:
    """Persist a finished simulation to ``obs_dir`` (mirrors the C++ output layout)."""
    obs_dir = Path(obs_dir)
    obs_dir.mkdir(parents=True, exist_ok=True)

    np.save(obs_dir / "internal_timestamps.npy", np.asarray(sim.solver.get_internal_time_stamps(), dtype=np.float64))

    npz_path = str(obs_dir / "unit_operations.npz")
    for name, obs in sim.observers.items():
        obs.save_to_npz(npz_path, name, "a")

    np.save(obs_dir / "pipe_outlet_middle_cell_pH_activity.npy", sim.outlet_pH(activity=True))
    np.save(obs_dir / "pipe_outlet_middle_cell_pH_concentration.npy", sim.outlet_pH(activity=False))
    return obs_dir

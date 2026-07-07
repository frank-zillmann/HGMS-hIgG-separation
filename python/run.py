"""Batch runner: solve the default-recipe experiment, save observations and figures to disk.

Command-line alternative to ``ui.py``. Run with ``python run.py``.
"""

from datetime import datetime
from pathlib import Path
from typing import Optional, Union

import fs3
from experiment import build_experiment
from plot_fractions import build_fractions_figure, default_datasets, load_reference_fractions
from plot_pH import build_ph_figure, load_experimental
from plot_time_step_sizes import build_time_step_figure
from observers import save_observations
from recipe import phase_transitions

# Where runs are saved (and where animate_unit_operations.py looks for one).
DEFAULT_RUN_DIR = Path(__file__).resolve().parent.parent / "data" / "runs"


def run(
    *,
    discretization_factor: float = 0.5,
    tau_reaction: float = 0.1,
    solver_type=fs3.SolverType.ERK,
    timeout_seconds: float = float("inf"),
    output_dir: Optional[Union[str, Path]] = None,
):
    experiment = build_experiment(
        discretization_factor=discretization_factor,
        tau_reaction=tau_reaction,
        solver_type=solver_type,
        timeout_seconds=timeout_seconds,
    )
    print(f"Total duration: {experiment.total_duration} s, {experiment.n_obs_time_steps} snapshots.")
    experiment.solve()
    print(f"Solved in {experiment.solver.get_solve_time():.1f} s ({len(experiment.solver.get_internal_time_stamps())} internal steps).")

    output_dir = Path(output_dir) if output_dir is not None else DEFAULT_RUN_DIR / datetime.now().strftime("run_%Y-%m-%d_%H-%M-%S")
    save_observations(experiment, output_dir / "obs")

    span = (0.0, experiment.total_duration)
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    build_ph_figure(
        experiment.times(), experiment.outlet_pH(True), experiment.outlet_pH(False),
        experimental=load_experimental(),
        transitions=phase_transitions(experiment.recipe), x_range=span,
    ).write_image(str(plots_dir / "plot_pH.pdf"))
    build_fractions_figure(
        default_datasets(experiment.fraction_masses(), load_reference_fractions())
    ).write_image(str(plots_dir / "plot_fractions.pdf"))
    build_time_step_figure(*experiment.step_sizes(), x_range=span).write_image(str(plots_dir / "plot_time_step_sizes.png"), scale=3)

    print(f"Saved observations and plots to {output_dir}")
    return experiment


if __name__ == "__main__":
    run()

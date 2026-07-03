"""Batch runner: solve an experiment, save observations and figures to disk.

This is the command-line alternative to ``ui.py``. Run with ``python run.py``.
"""

from datetime import datetime
from pathlib import Path
from typing import Optional, Sequence

from experiment import build_simulation
from figures import result_figures
from observers import save_observations
from plot_pH import load_experimental
from recipe import RecipeStep

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
RUNS_DIR = DATA_DIR / "runs"
EXPERIMENTAL_PATH = DATA_DIR / "experimental_data.ods"

# A previously saved run to inspect standalone (e.g. by animate_unit_operations.py).
DEFAULT_RUN_DIR = DATA_DIR / "run_fitted_no_MNP_hydroxyl_reactions_ERK_for_time_steps"


def run(
    recipe: Optional[Sequence[RecipeStep]] = None,
    *,
    discretization_factor: float = 0.5,
    run_tag: Optional[str] = None,
    output_base_dir: Optional[Path] = None,
    save: bool = True,
):
    sim = build_simulation(recipe, discretization_factor=discretization_factor)
    print(f"Total duration: {sim.total_duration} s, {sim.n_obs_time_steps} snapshots.")
    sim.solve()
    print(f"Solved in {sim.solver.get_solve_time():.1f} s ({len(sim.solver.get_internal_time_stamps())} internal steps).")

    results = sim.results()
    if save:
        run_dir = Path(output_base_dir or RUNS_DIR) / (run_tag or datetime.now().strftime("run_%Y-%m-%d_%H-%M-%S"))
        save_observations(sim, run_dir / "obs")

        plots_dir = run_dir / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        figures = result_figures(results, load_experimental(str(EXPERIMENTAL_PATH)))
        figures["pH"].write_image(str(plots_dir / "plot_pH.pdf"))
        figures["fractions"].write_image(str(plots_dir / "plot_fractions.pdf"))
        figures["time_steps"].write_image(str(plots_dir / "plot_time_step_sizes.png"), scale=3)
        print(f"Saved observations and plots to {run_dir}")

    return sim, results


if __name__ == "__main__":
    run()

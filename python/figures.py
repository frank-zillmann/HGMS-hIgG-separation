"""Turn a ``Results`` object into the standard set of Plotly figures (shared by run.py and ui.py)."""

from typing import Dict, Optional

import pandas as pd
import plotly.graph_objects as go

from plot_fractions import EXPERIMENT, build_fractions_figure
from plot_pH import build_ph_figure
from plot_time_step_sizes import build_time_step_figure


def result_figures(results, experimental: Optional[pd.DataFrame] = None) -> Dict[str, go.Figure]:
    x_range = (0.0, float(results.times[-1])) if len(results.times) else None
    return {
        "pH": build_ph_figure(
            results.times, results.ph_activity, results.ph_concentration,
            experimental=experimental, transitions=results.transitions, x_range=x_range,
        ),
        "fractions": build_fractions_figure({"Simulation": results.fraction_masses, "Experiment": EXPERIMENT}),
        "time_steps": build_time_step_figure(results.step_times, results.step_sizes),
    }

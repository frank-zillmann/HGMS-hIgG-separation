"""Fraction-mass figure builder (Plotly stacked bar).

The measured "Experiment" values (elutions only) from the paper figure can be added as an
extra scenario; Feed/Wash are unknown there and drawn as gaps.

Units: grams [g]
"""

from typing import Dict, Sequence

import numpy as np
import plotly.graph_objects as go

from recipe import FRACTIONS

# Elutions aligned first, then Feed and Wash stacked above.
DEFAULT_PLOT_ORDER = FRACTIONS[-5:] + FRACTIONS[:5]

# Light/bright colors so the value labels stay readable.
FRACTION_COLORS = {
    "Feed 1": "#d65f5f", "Feed 2": "#e07b7b",
    "Wash 1": "#5aa469", "Wash 2": "#7cbf7a", "Wash 3": "#a3d9a5",
    "Elution 1": "#6fa3d5", "Elution 2": "#85b5de", "Elution 3": "#9dc4e6",
    "Elution 4": "#b5d3ec", "Elution 5": "#cfe2f3",
}

EXPERIMENT = {
    "Elution 1": 0.32, "Elution 2": 0.63, "Elution 3": 0.42, "Elution 4": 0.32, "Elution 5": 0.06,
}

# Reference datasets from the paper/Matlab (elutions only; Feed/Wash unknown).
# Uncomment and merge into ``datasets`` to overlay them.
# REFERENCE_DATASETS = {
#     "Matlab, fitted Langmuir, main simulation":
#         {"Elution 1": 0.32, "Elution 2": 0.63, "Elution 3": 0.42, "Elution 4": 0.29, "Elution 5": 0.09},
#     "Matlab, isotherm Langmuir, main simulation":
#         {"Elution 1": 0.22, "Elution 2": 0.43, "Elution 3": 0.35, "Elution 4": 0.25, "Elution 5": 0.13},
#     "Matlab, fitted Langmuir, no slurry formation":
#         {"Elution 1": 0.75, "Elution 2": 0.43, "Elution 3": 0.30, "Elution 4": 0.20, "Elution 5": 0.09},
#     "Matlab, fitted Langmuir, no MNP hydroxyl reactions":
#         {"Elution 1": 0.36, "Elution 2": 0.68, "Elution 3": 0.54, "Elution 4": 0.21, "Elution 5": 0.06},
#     "Matlab, fitted Langmuir, no Washburn model":
#         {"Elution 1": 0.17, "Elution 2": 0.56, "Elution 3": 0.52, "Elution 4": 0.33, "Elution 5": 0.10},
# }


def build_fractions_figure(
    datasets: Dict[str, Dict[str, float]],
    plotted_fraction_names: Sequence[str] = DEFAULT_PLOT_ORDER,
    component: str = "hIgG",
) -> go.Figure:
    """One stacked bar per scenario; ``datasets`` maps scenario -> {fraction: mass_g}."""
    scenarios = list(datasets)
    fig = go.Figure()
    for fraction in plotted_fraction_names:
        values = [datasets[s].get(fraction, np.nan) for s in scenarios]
        if all(np.isnan(v) for v in values):
            continue
        fig.add_bar(
            name=fraction,
            x=scenarios,
            y=values,
            marker_color=FRACTION_COLORS.get(fraction),
            marker_line=dict(color="white", width=1.5),
            text=[f"{v:.2g}" if not np.isnan(v) else "" for v in values],
            textposition="inside",
            insidetextanchor="middle",
            textfont=dict(color="black", size=11),
        )
    fig.update_layout(
        barmode="stack",
        bargap=0.55,
        title=f"Fractions of {component}",
        template="plotly_white",
        yaxis_title=f"Mass {component} [g]",
        legend=dict(title="Fraction"),
        margin=dict(l=60, r=20, t=60, b=50),
    )
    return fig

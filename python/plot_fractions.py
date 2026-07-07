"""Fraction-mass figure builder (Plotly stacked bar).

The measured "Experiment" values (elutions only) from the paper figure can be added as an
extra scenario; Feed/Wash are unknown there and drawn as gaps.

Units: grams [g]
"""

from typing import Dict, Sequence

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from recipe import FRACTIONS

# Elutions aligned first, then Feed and Wash stacked above.
DEFAULT_PLOT_ORDER = FRACTIONS[-5:] + FRACTIONS[:5]

# Reference datasets (data/reference_fractions.csv) shown to the right of the simulation by default.
DEFAULT_REFERENCES = ["Experiment", "Matlab, isotherm Langmuir, main simulation"]

# Light/bright colors so the value labels stay readable.
FRACTION_COLORS = {
    "Feed 1": "#d65f5f", "Feed 2": "#e07b7b",
    "Wash 1": "#5aa469", "Wash 2": "#7cbf7a", "Wash 3": "#a3d9a5",
    "Elution 1": "#6fa3d5", "Elution 2": "#85b5de", "Elution 3": "#9dc4e6",
    "Elution 4": "#b5d3ec", "Elution 5": "#cfe2f3",
}


def load_reference_fractions(path: str) -> Dict[str, Dict[str, float]]:
    """Load {scenario: {fraction: mass_g}} from a CSV (scenario index, fraction columns)."""
    df = pd.read_csv(path, index_col="scenario")
    return {scenario: {frac: v for frac, v in row.items() if pd.notna(v)} for scenario, row in df.iterrows()}


def default_datasets(simulation_masses: Dict[str, float], references: Dict[str, Dict[str, float]]) -> Dict[str, Dict[str, float]]:
    """Simulation first, then the default reference scenarios (in order) that are available."""
    datasets = {"Simulation": simulation_masses}
    for name in DEFAULT_REFERENCES:
        if name in references:
            datasets[name] = references[name]
    return datasets


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

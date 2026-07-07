"""pH figure builder (Plotly)."""

from pathlib import Path
from typing import Optional, Sequence, Tuple

import pandas as pd
import plotly.graph_objects as go

EXPERIMENTAL_PATH = Path(__file__).resolve().parent.parent / "data" / "experimental_pH.csv"


def build_ph_figure(
    times: Sequence[float],
    ph_activity: Optional[Sequence[float]] = None,
    ph_concentration: Optional[Sequence[float]] = None,
    experimental: Optional[pd.DataFrame] = None,
    transitions: Optional[Sequence[float]] = None,
    x_range: Optional[Tuple[float, float]] = None,
    title: str = "pH in Outlet Pipe",
) -> go.Figure:
    """Build the pH figure. Accepts partial data so it can be redrawn live during a run."""
    fig = go.Figure()

    for t in transitions or []:
        fig.add_vline(x=t, line_width=0.5, line_color="gray", opacity=0.3)

    if experimental is not None:
        fig.add_scatter(
            x=experimental["time_s"], y=experimental["pH"] - 0.35,
            name="Experimental", mode="lines", line=dict(color="black", width=1.5),
        )
    if len(times):
        if ph_activity is not None and len(ph_activity):
            fig.add_scatter(
                x=times, y=ph_activity,
                name="Simulation (activity)", mode="lines", line=dict(color="#555555", width=1.5),
            )
        if ph_concentration is not None and len(ph_concentration):
            fig.add_scatter(
                x=times, y=ph_concentration,
                name="Simulation (concentration)", mode="lines",
                line=dict(color="#999999", width=1.5, dash="dash"),
            )

    fig.update_layout(
        title=title,
        template="plotly_white",
        xaxis=dict(title="Time [s]", dtick=500, range=list(x_range) if x_range else None),
        yaxis=dict(title="pH [-]", dtick=0.5),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        margin=dict(l=60, r=20, t=60, b=50),
    )
    return fig


def load_experimental() -> Optional[pd.DataFrame]:
    """Read the optional experimental pH reference (CSV with time_s,pH); None if unavailable."""
    if not EXPERIMENTAL_PATH.exists():
        print(f"Warning: Experimental data not found at {EXPERIMENTAL_PATH}.")
        return None
    return pd.read_csv(EXPERIMENTAL_PATH)[["time_s", "pH"]].sort_values("time_s")

import os
from typing import List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from model import path_to_obs as shared_path_to_obs
from model import path_to_save as shared_path_to_save


def build_ph_figure(
    times: Sequence[float],
    ph_activity: Sequence[float],
    ph_concentration: Optional[Sequence[float]] = None,
    experimental: Optional[pd.DataFrame] = None,
    transitions: Optional[Sequence[float]] = None,
    x_range: Optional[Tuple[float, float]] = None,
    title: str = "pH in Process Chamber",
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
    fig.add_scatter(
        x=times, y=ph_activity,
        name="Simulation (activity)", mode="lines", line=dict(color="#555555", width=1.5),
    )
    if ph_concentration is not None:
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


def load_experimental(
    experimental_path: str, t_from: float, t_until: float
) -> Optional[pd.DataFrame]:
    if not os.path.exists(experimental_path):
        print(f"Warning: Experimental data not found at {experimental_path}.")
        return None
    try:
        df = pd.read_excel(experimental_path, engine="odf")
    except Exception as exc:  # missing odf engine, bad file, ...
        print(f"Warning: Could not read experimental data ({exc}).")
        return None
    df = df[["time_s", "pH"]].sort_values("time_s")
    return df[(df["time_s"] >= t_from) & (df["time_s"] <= t_until)]


def plot_pH(
    path_to_obs: Optional[str] = None,
    path_to_save: Optional[str] = None,
    experimental_path: str = "data/experimental_data.ods",
    t_plot_from: float = 0.0,
    t_plot_until: float = 6600.0,
    phase_transitions: Optional[List[float]] = None,
    include_experimental: bool = True,
    save_plots: bool = True,
) -> Tuple[go.Figure, Optional[str]]:
    path_to_obs = path_to_obs or shared_path_to_obs
    path_to_save = path_to_save or shared_path_to_save

    activity = np.load(os.path.join(path_to_obs, "pipe_outlet_middle_cell_pH_activity.npy"))
    concentration = np.load(os.path.join(path_to_obs, "pipe_outlet_middle_cell_pH_concentration.npy"))
    times = np.arange(0.0, len(activity) * 2.0, 2.0)

    mask = (times >= t_plot_from) & (times <= t_plot_until)
    times, activity, concentration = times[mask], activity[mask], concentration[mask]

    experimental = None
    if include_experimental:
        experimental = load_experimental(experimental_path, t_plot_from, t_plot_until)

    transitions = [t for t in (phase_transitions or []) if t_plot_from <= t <= t_plot_until]

    fig = build_ph_figure(
        times, activity, concentration, experimental, transitions,
        x_range=(t_plot_from, t_plot_until),
    )

    outfile = None
    if save_plots:
        os.makedirs(path_to_save, exist_ok=True)
        outfile = os.path.join(path_to_save, "plot_pH.pdf")
        fig.write_image(outfile)

    return fig, outfile


if __name__ == "__main__":
    plot_pH()

import os
from typing import Optional, Sequence, Tuple

import numpy as np
import plotly.graph_objects as go

from model import path_to_obs as shared_path_to_obs
from model import path_to_save as shared_path_to_save


def build_time_step_figure(
    times: Sequence[float],
    step_sizes: Sequence[float],
    title: str = "Time Step Sizes taken by the Solver",
) -> go.Figure:
    # Scattergl (WebGL) keeps the tens of thousands of points responsive.
    fig = go.Figure(
        go.Scattergl(x=times, y=step_sizes, mode="markers", marker=dict(size=3, opacity=0.5))
    )
    fig.update_layout(
        title=title,
        template="plotly_white",
        xaxis_title="Time [s]",
        yaxis_title="Step size [s]",
        margin=dict(l=60, r=20, t=60, b=50),
    )
    return fig


def plot_time_step_sizes(
    path_to_obs: Optional[str] = None,
    path_to_save: Optional[str] = None,
    t_start: float = 0.0,
    t_end: float = 6600.0,
    save_plots: bool = True,
) -> Tuple[go.Figure, Optional[str]]:
    path_to_obs = path_to_obs or shared_path_to_obs
    path_to_save = path_to_save or shared_path_to_save

    timestamps = np.load(os.path.join(path_to_obs, "internal_timestamps.npy"))
    dt = np.diff(timestamps)
    timestamps = timestamps[:-1]

    mask = (timestamps >= t_start) & (timestamps <= t_end)
    fig = build_time_step_figure(timestamps[mask], dt[mask])

    outfile = None
    if save_plots:
        os.makedirs(path_to_save, exist_ok=True)
        outfile = os.path.join(path_to_save, "plot_time_step_sizes.png")
        fig.write_image(outfile, scale=3)

    return fig, outfile


if __name__ == "__main__":
    plot_time_step_sizes()

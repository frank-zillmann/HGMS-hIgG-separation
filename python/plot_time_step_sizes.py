"""Solver time-step-size figure builder (Plotly)."""

from typing import Optional, Sequence, Tuple

import plotly.graph_objects as go


def build_time_step_figure(
    times: Sequence[float],
    step_sizes: Sequence[float],
    x_range: Optional[Tuple[float, float]] = None,
    title: str = "Time Step Sizes taken by the Solver",
) -> go.Figure:
    # Scattergl (WebGL) keeps the tens of thousands of points responsive.
    fig = go.Figure(
        go.Scattergl(x=times, y=step_sizes, mode="markers", marker=dict(size=3, opacity=0.5))
    )
    fig.update_layout(
        title=title,
        template="plotly_white",
        xaxis=dict(title="Time [s]", range=list(x_range) if x_range else None),
        yaxis_title="Step size [s]",
        margin=dict(l=60, r=20, t=60, b=50),
    )
    return fig

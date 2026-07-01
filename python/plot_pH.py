import os
from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from model import path_to_obs as shared_path_to_obs
from model import path_to_save as shared_path_to_save


def plot_pH(
    path_to_obs: Optional[str] = None,
    path_to_save: Optional[str] = None,
    experimental_path: str = "data/experimental_data.ods",
    t_plot_from: float = 0.0,
    t_plot_until: float = 6600.0,
    phase_transitions: Optional[List[float]] = None,
    include_experimental: bool = True,
    save_plots: bool = True,
    close_fig: bool = True,
) -> Tuple[plt.Figure, Optional[str]]:
    if path_to_obs is None:
        path_to_obs = shared_path_to_obs
    if path_to_save is None:
        path_to_save = shared_path_to_save

    data_experimental = None
    if include_experimental and os.path.exists(experimental_path):
        df_experimental = pd.read_excel(experimental_path, engine="odf")
        data_experimental = df_experimental[["time_s", "pH"]].sort_values("time_s")
        data_experimental = data_experimental[
            (data_experimental["time_s"] >= t_plot_from)
            & (data_experimental["time_s"] <= t_plot_until)
        ]
    elif include_experimental:
        print(f"Warning: Experimental data not found at {experimental_path}.")

    pH_simulation_activity = np.load(os.path.join(path_to_obs, "pipe_outlet_middle_cell_pH_activity.npy"))
    pH_simulation_concentration = np.load(os.path.join(path_to_obs, "pipe_outlet_middle_cell_pH_concentration.npy"))
    time_stamps_simulation = np.arange(0.0, len(pH_simulation_activity) * 2.0, 2.0)
    assert len(time_stamps_simulation) == len(pH_simulation_activity) == len(pH_simulation_concentration), (
        "Simulation data length mismatch."
    )

    mask_simulation = (time_stamps_simulation >= t_plot_from) & (time_stamps_simulation <= t_plot_until)
    time_stamps_simulation = time_stamps_simulation[mask_simulation]
    pH_simulation_activity = pH_simulation_activity[mask_simulation]
    pH_simulation_concentration = pH_simulation_concentration[mask_simulation]

    fig = plt.figure(figsize=(20, 10))
    if data_experimental is not None:
        plt.plot(
            data_experimental["time_s"],
            data_experimental["pH"] - 0.35,
            linewidth=1.5,
            color="black",
            label="Experimental",
        )
    plt.plot(
        time_stamps_simulation,
        pH_simulation_activity,
        linewidth=1.5,
        color="grey",
        label="Simulation (activity based)",
    )
    plt.plot(
        time_stamps_simulation,
        pH_simulation_concentration,
        linewidth=1.5,
        color="grey",
        linestyle="--",
        label="Simulation (concentration based)",
    )

    if phase_transitions is None:
        phase_transitions = [
            0, 1199, 1233, 1329, 1362, 1412, 1477, 1528, 1597, 1793, 1799,
            1849, 1914, 1965, 2034, 2230, 2244, 2294, 2359, 2410, 2479, 2675,
            2688, 2740, 2830, 3130, 3190, 3390, 3397, 3449, 3539, 3839, 3899,
            4099, 4108, 4160, 4250, 4550, 4610, 4810, 4820, 4872, 4962, 5262,
            5322, 5522, 5676, 5728, 5818, 6118, 6178, 6378, 6548, 6600,
        ]

    phase_transitions_filtered = [t for t in phase_transitions if t_plot_from <= t <= t_plot_until]

    ax = plt.gca()
    ax.set_xlim(t_plot_from, t_plot_until)
    ax.autoscale(axis="y")
    y_min, y_max = ax.get_ylim()
    ax.set_yticks(np.arange(np.floor(y_min * 2) / 2, np.ceil(y_max * 2) / 2 + 0.5, 0.5))

    x_ticks_width = 250
    x_ticks = np.arange(np.ceil(t_plot_from / x_ticks_width) * x_ticks_width, t_plot_until + 1, x_ticks_width)
    ax.set_xticks(x_ticks)

    for transition in phase_transitions_filtered:
        ax.axvline(x=transition, color="gray", alpha=0.3, linewidth=0.5, linestyle="-")

    plt.grid(True, axis="y", alpha=0.3)

    plt.xlabel("Time [s]")
    plt.ylabel("pH [-]")
    plt.title("pH in Process Chamber")
    plt.legend()
    plt.tight_layout()

    outfile = None
    if save_plots:
        os.makedirs(path_to_save, exist_ok=True)
        outfile = os.path.join(path_to_save, "plot_pH.pdf")
        plt.savefig(outfile, bbox_inches="tight")

    if close_fig:
        plt.close(fig)

    return fig, outfile


if __name__ == "__main__":
    plot_pH()


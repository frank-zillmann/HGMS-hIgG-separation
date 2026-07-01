"""
Plot script for fractions.

Adds a convenient, hard-coded dataset with the exact numbers from the
paper figure (Experiment, Simulation 1, Simulation 2, Simulation 1A, 1B, 1C).
Feed and Wash values are not available in the figure and are therefore set to
np.nan so you can access a uniform structure without special cases.

Units: grams [g]
"""

import os
from typing import Dict, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

from model import build_component_system
from model import path_to_obs_dict as shared_path_to_obs_dict
from model import path_to_save as shared_path_to_save

# ---- Static values copied from the figure (all in grams) --------------------
# Fraction ordering keys used throughout (Feed and Wash intentionally NaN)
all_fraction_names = [
    "Feed 1",
    "Feed 2",
    "Wash 1",
    "Wash 2",
    "Wash 3",
    "Elution 1",
    "Elution 2",
    "Elution 3",
    "Elution 4",
    "Elution 5"
]

# Only plot elution fractions in stacked bars
#plotted_fraction_names = all_fraction_names[-5:]

# Plot all fractions including Feed and Wash
# plotted_fraction_names = all_fraction_names

# Plot Elutions aligned and Feed and Wash above
plotted_fraction_names = all_fraction_names[-5:] + all_fraction_names[:5]

# Consistent, light/bright colors so labels remain readable
# - Feed/Wash use soft pastels
# - Elutions use a light grayscale gradient from slightly darker (E1) to very light (E5)
fraction_colors = {
    # Feed (muted reds)
    "Feed 1": {"face": "#d65f5f", "edge": "#000000"},
    "Feed 2": {"face": "#e07b7b", "edge": "#000000"},

    # Wash (muted greens)
    "Wash 1": {"face": "#5aa469", "edge": "#888888"},
    "Wash 2": {"face": "#7cbf7a", "edge": "#888888"},
    "Wash 3": {"face": "#a3d9a5", "edge": "#888888"},

    # Elution (blue gradient: E1 moderately saturated -> E5 light)
    "Elution 1": {"face": "#6fa3d5", "edge": None},
    "Elution 2": {"face": "#85b5de", "edge": None},
    "Elution 3": {"face": "#9dc4e6", "edge": None},
    "Elution 4": {"face": "#b5d3ec", "edge": None},
    "Elution 5": {"face": "#cfe2f3", "edge": None},
}

def plot_fractions(
    path_to_obs_dict: Optional[Dict[str, str]] = None,
    path_to_save: Optional[str] = None,
    component: str = "hIgG",
    plotted_fraction_names: Optional[list] = None,
    save_plots: bool = True,
    close_fig: bool = True,
) -> Tuple[plt.Figure, Optional[str]]:
    if path_to_obs_dict is None:
        path_to_obs_dict = shared_path_to_obs_dict
    if path_to_save is None:
        path_to_save = shared_path_to_save

    factor_raw_to_g = 1000
    component_idx = build_component_system().get_idx(component)
    time_idx = -1  # Last time point

    if plotted_fraction_names is None:
        plotted_fraction_names = all_fraction_names[-5:] + all_fraction_names[:5]

    all_datasets: Dict[str, Dict[str, float]] = {}

    for sim_name, obs_path in path_to_obs_dict.items():
        npz_path = os.path.join(obs_path, "unit_operations.npz")
        if not os.path.exists(npz_path):
            print(f"Warning: Missing simulation data at {npz_path}; skipping '{sim_name}'.")
            continue
        unit_operation_data = np.load(npz_path)

        values_g = np.array([
            unit_operation_data[name][time_idx, 0, component_idx] for name in all_fraction_names
        ]) * factor_raw_to_g

        all_datasets[sim_name] = {name: val for name, val in zip(all_fraction_names, values_g)}

    _nan_block = {
        "Feed 1": np.nan,
        "Feed 2": np.nan,
        "Wash 1": np.nan,
        "Wash 2": np.nan,
        "Wash 3": np.nan,
    }

    all_datasets["Experiment"] = {
        **_nan_block,
        "Elution 1": 0.32,
        "Elution 2": 0.63,
        "Elution 3": 0.42,
        "Elution 4": 0.32,
        "Elution 5": 0.06,
    }
# all_datasets["Matlab, fitted Langmuir, main simulation"] = {
#         **_nan_block,
#         "Elution 1": 0.32,
#         "Elution 2": 0.63,
#         "Elution 3": 0.42,
#         "Elution 4": 0.29,
#         "Elution 5": 0.09,
#     }
# all_datasets["Matlab, isotherm Langmuir, main simulation"] = {
#         **_nan_block,
#         "Elution 1": 0.22,
#         "Elution 2": 0.43,
#         "Elution 3": 0.35,
#         "Elution 4": 0.25,
#         "Elution 5": 0.13,
#     }
# all_datasets["Matlabm, fitted Langmuir, no slurry formation"] = {
#         **_nan_block,
#         "Elution 1": 0.75,
#         "Elution 2": 0.43,
#         "Elution 3": 0.30,
#         "Elution 4": 0.20,
#         "Elution 5": 0.09,
#     }
# all_datasets["Matlab, fitted Langmuir, no MNP hydroxyl reactions"] = {
#         **_nan_block,
#         "Elution 1": 0.36,
#         "Elution 2": 0.68,
#         "Elution 3": 0.54,
#         "Elution 4": 0.21,
#         "Elution 5": 0.06,
#     }
# all_datasets["Matlab, fitted Langmuir, no Washburn model"] = {
#         **_nan_block,
#         "Elution 1": 0.17,
#         "Elution 2": 0.56,
#         "Elution 3": 0.52,
#         "Elution 4": 0.33,
#         "Elution 5": 0.10,
#     }

    print(f"All datasets: {all_datasets}")

    fig_width = max(8, int(1.2 * len(all_datasets) + 4))
    fig, ax = plt.subplots(figsize=(fig_width, 6))

    spacing_factor = 0.5
    x = np.arange(len(all_datasets)) * spacing_factor
    bar_width = 0.4

    for si, scenario in enumerate(all_datasets.keys()):
        bottom = 0.0
        for fraction_name in plotted_fraction_names:
            value = all_datasets[scenario][fraction_name]
            if np.isnan(value):
                continue
            color_info = fraction_colors.get(fraction_name, None)
            if isinstance(color_info, dict):
                face_color = color_info["face"]
                edge_color = color_info["edge"] if color_info["edge"] is not None else "#ffffff"
            else:
                face_color = color_info
                edge_color = "#ffffff"
            ax.bar(
                x[si],
                value,
                width=bar_width,
                bottom=bottom,
                color=face_color,
                edgecolor=edge_color,
                linewidth=1.5,
                label=fraction_name if si == 0 else None,
            )
            y = bottom + value / 2.0
            ax.text(x[si], y, f"{value:.2g}", ha="center", va="center", color="black", fontsize=9)
            print(f"{scenario} — {fraction_name}: {value} g")
            bottom += value

    ax.set_xticks(x)
    ax.set_xticklabels(
        list(all_datasets.keys()),
        rotation=-10,
        ha="left",
        rotation_mode="anchor",
    )
    ax.tick_params(axis="x", pad=4)
    ax.margins(x=0.0)

    left_pad = 0.1
    right_pad = 0.6
    ax.set_xlim(x[0] - bar_width / 2 - left_pad, x[-1] + bar_width / 2 + right_pad)
    ax.set_ylabel(f"Mass {component} [g]")
    ax.set_title(f"Fractions of {component}")
    ax.legend(title="Fraction", loc="upper right")
    plt.tight_layout()

    outfile = None
    if save_plots:
        os.makedirs(path_to_save, exist_ok=True)
        outfile = os.path.join(path_to_save, "plot_fractions.pdf")
        fig.savefig(outfile, bbox_inches="tight")
        print(f"Saved to {outfile}!")

    if close_fig:
        plt.close(fig)

    return fig, outfile


if __name__ == "__main__":
    plot_fractions()

import os
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

from shared_data import path_to_obs as shared_path_to_obs
from shared_data import path_to_save as shared_path_to_save
from shared_data import show_plots as shared_show_plots


def plot_time_step_sizes(
	path_to_obs: Optional[str] = None,
	path_to_save: Optional[str] = None,
	show_plots: Optional[bool] = None,
	t_start: float = 0.0,
	t_end: float = 6600.0,
	save_plots: bool = True,
	close_fig: bool = True,
) -> Tuple[plt.Figure, Optional[str]]:
	if path_to_obs is None:
		path_to_obs = shared_path_to_obs
	if path_to_save is None:
		path_to_save = shared_path_to_save
	if show_plots is None:
		show_plots = shared_show_plots

	timestamps = np.load(os.path.join(path_to_obs, "internal_timestamps.npy"))

	dt = np.diff(timestamps)
	timestamps = timestamps[0:-1]

	idx_start = np.abs(timestamps - t_start).argmin()
	idx_end = np.abs(timestamps - t_end).argmin()

	fig = plt.figure(figsize=(10, 6))
	plt.plot(timestamps[idx_start:idx_end], dt[idx_start:idx_end], marker=".", markersize=3, linestyle="")
	plt.xlabel("Time [s]")
	plt.ylabel("Step size [s]")
	plt.title("Time Step Sizes taken by the Solver")
	plt.grid(True)

	outfile = None
	if save_plots:
		os.makedirs(path_to_save, exist_ok=True)
		outfile = os.path.join(path_to_save, "plot_time_step_sizes.png")
		plt.savefig(outfile, bbox_inches="tight", dpi=300)

	if show_plots:
		plt.show()
	elif close_fig:
		plt.close(fig)

	return fig, outfile


if __name__ == "__main__":
	plot_time_step_sizes()

import numpy as np
import matplotlib.pyplot as plt
import os
from shared_data import path_to_obs, path_to_save, show_plots

# Load internal time stamps
timestamps =  np.load(path_to_obs + '/internal_timestamps.npy')

# Compute step sizes
dt = np.diff(timestamps)
timestamps = timestamps[0:-1]  # Align timestamps with dt

# Set time range to plot
#range = (0, -1)  # All

t_start = 0
t_end = 6600
range = (np.abs(timestamps - t_start).argmin(), np.abs(timestamps - t_end).argmin())

plt.figure(figsize=(10, 6))
plt.plot(timestamps[range[0]:range[-1]], dt[range[0]:range[-1]], marker='.', markersize=3, linestyle="")
plt.xlabel('Time [s]')
plt.ylabel('Step size [s]')
plt.title('Time Step Sizes taken by the Solver')
plt.grid(True)

# Save to configured results directory as PNG (not as PDF to avoid too large files)
os.makedirs(path_to_save, exist_ok=True)
outfile = os.path.join(path_to_save, "plot_time_step_sizes.png")
plt.savefig(outfile, bbox_inches="tight", dpi = 300)

# Only show plots if configured
if show_plots:
	plt.show()
else:
	plt.close('all')

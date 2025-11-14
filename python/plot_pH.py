import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from shared_data import idx, path_to_obs, time_stamps, path_to_save, show_plots

t_plot_from = 0  # seconds
t_plot_until = 6600  # seconds

# Minimal read: assumes well-formatted ODS with columns: time_s, pH
df_experimental = pd.read_excel("data/experimental_data.ods", engine="odf")  # first sheet by default
# Keep only the columns we need, drop missing rows, and sort by time
data_experimental = df_experimental[["time_s", "pH"]].sort_values("time_s")
# Filter experimental data by time range
data_experimental = data_experimental[(data_experimental["time_s"] >= t_plot_from) & 
                                       (data_experimental["time_s"] <= t_plot_until)]

# data = np.load(path_to_obs + "/some_snapshot.npy")
time_stamps_simulation = time_stamps
pH_simulation_activity = np.load(path_to_obs + "/pipe_outlet_middle_cell_pH_activity.npy")
pH_simulation_concentration = np.load(path_to_obs + "/pipe_outlet_middle_cell_pH_concentration.npy")
assert len(time_stamps_simulation) == len(pH_simulation_activity) == len(pH_simulation_concentration), \
    "Simulation data length mismatch."

# Filter simulation data by time range
mask_simulation = (time_stamps_simulation >= t_plot_from) & (time_stamps_simulation <= t_plot_until)
time_stamps_simulation = time_stamps_simulation[mask_simulation]
pH_simulation_activity = pH_simulation_activity[mask_simulation]
pH_simulation_concentration = pH_simulation_concentration[mask_simulation]

print(f"Activity based pH vs Concentration based pH at t = {time_stamps_simulation[500]}: {pH_simulation_activity[500]} vs {pH_simulation_concentration[500]}")

plt.figure(figsize=(20, 10))
# subtract -0.35 from experimental data as reported in Marit Schneiders Master Thesis
plt.plot(data_experimental["time_s"], data_experimental["pH"] - 0.35, linewidth=1.5, color="black", label="Experimental")
plt.plot(time_stamps_simulation, pH_simulation_activity, linewidth=1.5, color="grey", label="Simulation (activity based)")
plt.plot(time_stamps_simulation, pH_simulation_concentration, linewidth=1.5, color="grey", linestyle="--", label="Simulation (concentration based)")

# Add grid with custom x positions at phase transitions and 0.5 spacing on y-axis
phase_transitions = [0, 1199, 1233, 1329, 1362, 1412, 1477, 1528, 1597, 1793, 1799, 
                     1849, 1914, 1965, 2034, 2230, 2244, 2294, 2359, 2410, 2479, 2675, 
                     2688, 2740, 2830, 3130, 3190, 3390, 3397, 3449, 3539, 3839, 3899, 
                     4099, 4108, 4160, 4250, 4550, 4610, 4810, 4820, 4872, 4962, 5262, 
                     5322, 5522, 5676, 5728, 5818, 6118, 6178, 6378, 6548, 6600]

# Filter phase transitions to only those within the plot range
phase_transitions_filtered = [t for t in phase_transitions if t_plot_from <= t <= t_plot_until]

# Set major ticks with labels at evenly spaced intervals
ax = plt.gca()
ax.set_xlim(t_plot_from, t_plot_until)

# Let matplotlib determine the y-axis range, then set ticks every 0.5
ax.autoscale(axis='y')
y_min, y_max = ax.get_ylim()
ax.set_yticks(np.arange(np.floor(y_min * 2) / 2, np.ceil(y_max * 2) / 2 + 0.5, 0.5))

# Set x-axis ticks every x_ticks_width
x_ticks_width = 250
x_ticks = np.arange(np.ceil(t_plot_from / x_ticks_width) * x_ticks_width, t_plot_until + 1, x_ticks_width)
ax.set_xticks(x_ticks)

# Draw grid lines at phase transitions by using vlines
for transition in phase_transitions_filtered:
    ax.axvline(x=transition, color='gray', alpha=0.3, linewidth=0.5, linestyle='-')

# Draw horizontal grid lines
plt.grid(True, axis='y', alpha=0.3)

plt.xlabel("Time [s]")
plt.ylabel("pH [-]")
plt.title("pH in Process Chamber")
plt.legend()
plt.tight_layout()

# Save to configured results directory as PDF
os.makedirs(path_to_save, exist_ok=True)
outfile = os.path.join(path_to_save, "plot_pH.pdf")
plt.savefig(outfile, bbox_inches="tight")

# Only show plots if configured
if show_plots:
    plt.show()
else:
    plt.close('all')


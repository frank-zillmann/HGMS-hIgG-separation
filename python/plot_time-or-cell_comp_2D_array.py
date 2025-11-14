import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import os
from shared_data import idx, path_to_obs, path_to_save, show_plots

# data = np.load(path_to_obs + "/some_snapshot.npy")
data = np.load(path_to_obs + "/pc_outlet.npy")

print(f"Shape: {data.shape}")
data = np.squeeze(data)
print(f"Squeezed Shape: {data.shape}")
print(f"Data type: {data.dtype}")

# List all components to plot (excluding H2O if desired)
#components = [c for c in idx.keys() if c != "H₂O"]
components = ["H⁺", "MNP-OH", "MNP1", "MNP10", "hIgG", "MNP-hIgG"]

num_components = len(components)
fig, axes = plt.subplots(num_components, 1, figsize=(10, 2*num_components), sharex=True)

if num_components == 1:
    axes = [axes]

for ax, component in zip(axes, components):
    ax.plot(data[:, idx[component]], label=component)
    ax.set_ylabel(component)
    ax.legend(loc='upper right')
    print(f"Plotting component: {component} with index {idx[component]}")
    print(f"Data for {component}\n: {data[:, idx[component]]}")

axes[-1].set_xlabel('Cell index or Time index')
fig.suptitle('Mass per cell vectors (each component in own subplot)')
plt.tight_layout(rect=[0, 0.03, 1, 0.97])

# Save to configured results directory as PDF
os.makedirs(path_to_save, exist_ok=True)
outfile = os.path.join(path_to_save, "plot_time-or-cell_comp_2D_array.pdf")
plt.savefig(outfile, bbox_inches="tight")

# Only show plots if configured
if show_plots:
    plt.show()
else:
    plt.close('all')
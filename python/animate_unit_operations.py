"""
Animation script for Unit Operations layout.
- Accepts 6 numpy arrays (Pipe1, Loop, PC, PCsl, Dead, Pipe)
- Animates them in a layout mimicking the process diagram
- Option to plot each component in separate subplots or all together in one plot per Unit Operation
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import matplotlib.gridspec as gridspec
from shared_data import idx, path_to_obs, path_to_save, show_plots
import matplotlib.cm as cm


# --- CONFIGURATION ---
# Driven by shared_data: show when show_plots is True; otherwise, save.
show_animation = bool(show_plots)
save_animation = not show_animation
fps = 25

# Build output path under configured save directory
os.makedirs(path_to_save, exist_ok=True)
video_filename = os.path.join(path_to_save, "animation_unit_operations.mp4")

# Option: Plot each component in separate subplots (True) or all together (False)
separate_components = True

# --- LOAD DATA ---
unit_operation_data = np.load(path_to_obs + "/unit_operations.npz")

pipe_inlet_data = unit_operation_data["pipe_inlet"] # shape: (time, cells, components)
pipe_loop_data = unit_operation_data["pipe_loop"]
pc_liquid_data = unit_operation_data["pc_liquid"]
pc_slurry_data = unit_operation_data["pc_slurry"]
pipe_outlet_data = unit_operation_data["pipe_outlet"]

unit_ops = {
    "Inlet Pipe": [pipe_inlet_data, [], [], []],
    "Process Chamber Liquid": [pc_liquid_data, [], [], []],
    "Process Chamber Slurry": [pc_slurry_data, [], [], []],
    "Outlet Pipe": [pipe_outlet_data, [], [], []],
    "Loop Pipe": [pipe_loop_data, [], [], []]
 } # name : data, axes, lines, y_lims

# Choose components to plot
component_groups = [["H⁺", "OH⁻", "Na⁺", "Cl⁻", "TrisH⁺", "Tris", "AcH", "Ac⁻", "GlyH₂⁺", "GlyH", "Gly⁻"], 
                    ["MNP1", "MNP2", "MNP3", "MNP4", "MNP5", "MNP6", "MNP7", "MNP8", "MNP9", "MNP10"], # Without: "MNP-OH₂⁺", "MNP-OH", "MNP-O⁻", 
                    ["hIgG", "MNP-hIgG"]]
n_groups = len(component_groups)

color_groups = [
    # First group: different blue tones, H+ is darkest, MNP_OH2+/MNP_OH/MNP_O- are red/blue
    [
        "#060044",  # H+ (darkest blue)
        "#4335FF",   # OH-
        "#5B00A0",   # Na+
        "#AB3DFF",   # Cl-
        "#008055", # TrisH+
        "#00FFA6", # Tris
        "#498000",    # AcH
        "#91FF00",    # Ac-
        "#807D00",    # GlyH2+
        "#C9C500",    # GlyH
        "#F9F651"    # Gly⁻
    ],
    # [   
    #     "#6A0067",       # MNP_OH2+
    #     "#DA00D2",       # MNP_OH
    #     "#FF57F9"        # MNP_O-
    # ] 
    # + 
    ["#FCAFAF", "#FC8484", "#FB6464", "#FF4B4B", "#FF2929", "#C50000", "#930000", "#570000", "#360000", "#000000"],     # MNP1-MNP10, different red tones
    # Third group: hIgG orange, MNP_hIgG green
    ["#FF5500", "#0BE400"]
]

# --- FIGURE LAYOUT ---
fig = plt.figure(figsize=(15, 15))
gs = gridspec.GridSpec(3*n_groups, 3, height_ratios=n_groups*[1, 1, 1], width_ratios=[1, 2, 1])

# Axes
for i in range(len(component_groups)):
    unit_ops["Inlet Pipe"][1].append(fig.add_subplot(gs[1 * n_groups + i, 0]))
    unit_ops["Process Chamber Liquid"][1].append(fig.add_subplot(gs[1 * n_groups + i, 1]))
    unit_ops["Process Chamber Slurry"][1].append(fig.add_subplot(gs[2 * n_groups + i, 1]))
    unit_ops["Outlet Pipe"][1].append(fig.add_subplot(gs[1 * n_groups + i, 2]))
    unit_ops["Loop Pipe"][1].append(fig.add_subplot(gs[0 * n_groups + i, 1]))

for name, [data, group_ax, group_lines, group_lims] in unit_ops.items():
    for ax in group_ax:
        ax.set_ylabel('Mass per cell [kg]')
        component_names = ", ".join(component_groups[group_ax.index(ax)])
        if len(component_names) > 35:
            component_names = component_names[:32] + "..."
        ax.set_title(name + ": " + component_names)

# Prepare lines and y lims for each subplot
for name, [data, group_ax, group_lines, group_lims] in unit_ops.items():

    for i, components in enumerate(component_groups):
        component_lines = []
        group_y_min = float('inf')
        group_y_max = float('-inf')

        for j, comp in enumerate(components):
            color = color_groups[i][j]
            # Create line with dummy data so legend works
            line, = group_ax[i].plot([0], [0], label=comp, color=color)
            component_lines.append(line)

            comp_data = data[:, :, idx[comp]].flatten()
            y_min = np.min(comp_data)
            y_max = np.max(comp_data)
            # y_min = np.percentile(comp_data, 0)
            # y_max = np.percentile(comp_data, 100)
            if y_min < group_y_min:
                group_y_min = y_min
            if y_max > group_y_max:
                group_y_max = y_max

        group_lines.append(component_lines)
        group_lims.append((group_y_min, group_y_max))

fig.suptitle('Unit Operations Animation')
plt.tight_layout(rect=[0, 0.01, 1, 0.99])  # Make plots more tightly packed
plt.subplots_adjust(hspace=0.35, wspace=0.15)  # Increase vertical space between subplots

# Flatten all lines for animation
all_lines = []
for name, [data, group_ax, group_lines, group_lims] in unit_ops.items():
    for component_lines in group_lines:
        for line in component_lines:
            all_lines.append(line)

n_components = sum(len(g) for g in component_groups)
fig.legend(all_lines[0:n_components], [line.get_label() for line in all_lines[0:n_components]], loc='upper right', fontsize='large')

# --- ANIMATION FUNCTIONS ---
time_text = fig.text(0.98, 0.02, '', ha='right', va='bottom', color='black', fontsize=20)

def init():
    for line in all_lines:
        line.set_data([], [])
    time_text.set_text('')
    return all_lines + [time_text]

def animate(frame):
    for name, (data, group_ax, group_lines, group_lims) in unit_ops.items():
        x = np.arange(data.shape[1])
        for j, components in enumerate(component_groups):
            for k, comp in enumerate(components):
                y = data[frame, :, idx[comp]]
                group_lines[j][k].set_data(x, y)

            # Mirror x-axis for Loop axes only
            if name == "Loop Pipe":
                group_ax[j].set_xlim(data.shape[1]-1, 0)
            else:
                group_ax[j].set_xlim(0, data.shape[1]-1)
            group_ax[j].set_ylim(group_lims[j][0], group_lims[j][1])

    time_text.set_text(f'Time: {frame}/{timesteps}')
    return all_lines + [time_text]

# --- ANIMATION ---
# Use min time steps across all arrays
timesteps = unit_ops["Inlet Pipe"][0].shape[0]  # Assuming all have same time steps
ani = animation.FuncAnimation(
    fig, animate, frames=range(timesteps), init_func=init,
    blit=False, interval=int(1000 / fps), repeat=False
)

if show_animation:
    plt.show()

# --- SAVE OR SHOW ---
if save_animation:
    os.makedirs(path_to_save, exist_ok=True)
    if 'ffmpeg' in animation.writers:
        print("Using ffmpeg writer to save animation.")
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=fps, metadata=dict(artist='Me'), bitrate=1800)
        ani.save(video_filename, writer=writer)
        print(f"Animation saved to {video_filename}")
    elif 'pillow' in animation.writers:
        print("ffmpeg not found, using pillow writer to save animation as GIF.")
        gif_filename = os.path.splitext(video_filename)[0] + '.gif'
        Writer = animation.writers['pillow']
        writer = Writer(fps=fps)
        ani.save(gif_filename, writer=writer)
        print(f"Animation saved to {gif_filename} (GIF, not MP4)")
    else:
        print("No supported video writer (ffmpeg or pillow) found. Cannot save animation.")

# --- OPTIONAL: Draw arrows/lines to visually connect subplots (not animated) ---
# This can be done with fig.annotate or fig.lines if desired for static layout
# ...existing code...

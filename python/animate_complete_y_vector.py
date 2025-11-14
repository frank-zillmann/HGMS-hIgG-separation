import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import matplotlib.animation as animation
import os

from shared_data import idx, path_to_obs, path_to_save, show_plots

# Drive from shared_data: show if show_plots; otherwise save
show_animation = bool(show_plots)
save_animation = not show_animation
fps = 20 # Frames per second

# Output video filename under configured results directory
os.makedirs(path_to_save, exist_ok=True)
video_filename = os.path.join(path_to_save, "animation_complete_y_vector.mp4")

# Shape of data: (time_steps, cells, components)
data = np.load(path_to_obs + "/complete_y_vector.npy")

print(f"Shape: {data.shape}")
print(f"Data type: {data.dtype}")
print(f"Data: \n{data}")

# List all components to plot (excluding H2O if desired)
#components = [c for c in idx.keys() if c != "H₂O"]
#components = ["H⁺", "OH⁻", "Na⁺", "Cl⁻", "TrisH⁺", "Tris", "MNP1", "MNP2", "MNP3", "MNP4", "MNP5", "MNP6", "MNP7", "MNP8", "MNP9", "MNP10", "hIgG"]
components = ["H⁺", "Na⁺", "Cl⁻", "MNP1", "MNP10", "hIgG"]


# Animation with subplots for each component
num_components = len(components)
fig, axes = plt.subplots(num_components, 1, figsize=(10, 2*num_components), sharex=True)
if num_components == 1:
    axes = [axes]

# Adjustable frame rate (frames per second)
interval = int(1000 / fps)  # milliseconds per frame

lines = []
for ax, component in zip(axes, components):
    line, = ax.plot([], [], label=component)
    ax.set_ylabel(component)
    ax.legend(loc='upper right')
    lines.append(line)

axes[-1].set_xlabel('Cell index')
fig.suptitle('Concentration vectors (each component in own subplot)')
plt.tight_layout(rect=[0, 0.03, 1, 0.97])

# Add time index, start, and end time annotation in bottom right corner
start_time = 0
#end_time = 1800
end_time = data.shape[0] - 50

time_text = fig.text(0.98, 0.02, '', ha='right', va='bottom', fontsize=14, color='red')

def init():
    for line in lines:
        line.set_data([], [])
    time_text.set_text(f'Start: {start_time}   End: {end_time}')
    return lines + [time_text]

# Precompute y-limits for each component in the selected time range
y_lims = []
for component in components:
    comp_data = data[start_time:end_time+1, :, idx[component]]
    y_min = np.min(comp_data)
    y_max = np.max(comp_data)
    y_lims.append((y_min, y_max))

def animate(frame):
    for i, component in enumerate(components):
        # x: cell index, y: concentration at this time step
        x = np.arange(data.shape[1])
        y = data[frame, :, idx[component]]
        lines[i].set_data(x, y)
        axes[i].set_xlim(0, data.shape[1]-1)
        axes[i].set_ylim(y_lims[i][0], y_lims[i][1])
    time_text.set_text(f'Start: {start_time}   End: {end_time}   Time index: {frame}')
    return lines + [time_text]

# Use range(start_time, end_time+1) for frames
ani = animation.FuncAnimation(
    fig, animate, frames=range(start_time, end_time+1), init_func=init,
    blit=False, interval=interval, repeat=False
)

# Save animation to video if requested
if save_animation:
    os.makedirs(path_to_save, exist_ok=True)
    # Use ffmpeg writer if available, else fallback to pillow (gif)
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

# Show animation if requested
if show_animation:
    plt.show()

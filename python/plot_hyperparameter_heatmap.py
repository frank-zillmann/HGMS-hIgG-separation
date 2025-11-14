#!/usr/bin/env python3
"""
Script to visualize hyperparameter grid search results from the HGMS hIgG separation experiment.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from shared_data import path_to_save, show_plots

def load_and_plot_results(npz_file, value_type):
    """Load and plot a single heatmap of the selected value type."""
    # Load the data
    data = np.load(npz_file)
    print(f"All available file names: {data.files}")
    max_errors = data['max_errors']
    solve_times = data['solve_times']
    num_steps = data.get('num_steps')
    discretization_factors = data['discretization_factors']
    tau_values = data['tau_reaction_values']

    value_map = {
        'max_errors': (max_errors, 'Maximum Error', 'viridis', 'Max Error'),
        'solve_times': (solve_times, 'Solve Time (seconds)', 'plasma', 'Time (s)'),
        'num_steps': (num_steps, 'Number of Internal Steps', 'cividis', '# steps')
    }
    if value_type not in value_map or value_map[value_type][0] is None:
        print(f"Value type '{value_type}' not available in data.")
        return
    values, title, cmap, colorbar_label = value_map[value_type]

    # Plot single heatmap
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(values, aspect='auto', cmap=cmap, origin='lower')
    ax.set_xlabel(r'$\tau$ (reaction time constant)')
    ax.set_ylabel(r'discretization factor')
    ax.set_title(title)
    ax.set_xticks(range(len(tau_values)))
    ax.set_xticklabels([f'{v:.0e}' for v in tau_values], rotation=45)
    ax.set_yticks(range(len(discretization_factors)))
    ax.set_yticklabels([f'{v:.1f}' for v in discretization_factors])
    plt.colorbar(im, ax=ax, label=colorbar_label)

    # Add value text to each cell
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            val = values[i, j]
            if np.isnan(val):
                display_val = 'NaN'
            else:
                # Format based on value type
                if value_type == 'max_errors':
                    display_val = f'{val:.2e}'
                elif value_type == 'solve_times':
                    display_val = f'{val:.2f}'
                else:
                    display_val = f'{val:.0f}'
            ax.text(j, i, display_val, ha='center', va='center', color='white', fontsize=8, fontweight='bold')

    plt.tight_layout()
    os.makedirs(path_to_save, exist_ok=True)
    output_file = os.path.join(path_to_save, os.path.basename(npz_file).replace('.npz', f'_{value_type}_heatmap.pdf'))
    plt.savefig(output_file, bbox_inches='tight')
    print(f"Saved heatmap to: {output_file}")
    if show_plots:
        plt.show()
    else:
        plt.close('all')


if __name__ == "__main__":
    # Choose which value to plot: 'max_errors', 'solve_times', or 'num_steps'
    load_and_plot_results("data/hyperparameter_search_ADAMS.npz", "solve_times")
    load_and_plot_results("data/hyperparameter_search_ERK.npz", "solve_times")
    load_and_plot_results("data/hyperparameter_search_ADAMS.npz", "num_steps")
    load_and_plot_results("data/hyperparameter_search_ERK.npz", "num_steps")



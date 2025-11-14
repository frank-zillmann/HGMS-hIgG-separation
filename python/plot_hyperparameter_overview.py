#!/usr/bin/env python3
"""
Script to visualize hyperparameter grid search results from the HGMS hIgG separation experiment.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from shared_data import path_to_save, show_plots

def load_and_plot_results(npz_file):
    """Load and plot the hyperparameter search results."""
    
    # Load the data
    data = np.load(npz_file)

    print(f"All available file names: {data.files}")  # For debugging: print available arrays in the .npz file
    
    max_errors = data['max_errors']  # Shape: (n_disc, n_tau)
    solve_times = data['solve_times']  # Shape: (n_disc, n_tau)
    num_steps = data.get('num_steps')  # Shape: (n_disc, n_tau) if present
    discretization_factors = data['discretization_factors']  # Shape: (n_disc,)
    tau_values = data['tau_reaction_values']  # Shape: (n_tau,)
    
    print(f"Loaded data from: {npz_file}")
    print(f"Max errors (shape {max_errors.shape}):\n{max_errors}")
    print(f"Solve times (shape {solve_times.shape}):\n{solve_times}")
    print(f"discretization factors ({len(discretization_factors)}): {discretization_factors}")
    print(f"tau_reaction values ({len(tau_values)}): {tau_values}")
    
    # Create figure with subplots
    # Choose layout depending on whether num_steps is available
    has_steps = num_steps is not None
    fig_cols = 3 if has_steps else 2
    fig, axes = plt.subplots(2, fig_cols, figsize=(7*fig_cols, 12))
    
    # Plot 1: Max error heatmap
    ax = axes[0, 0]
    im1 = ax.imshow(max_errors, aspect='auto', cmap='viridis', origin='lower')
    ax.set_xlabel('tau index')
    ax.set_ylabel('discretization index')
    ax.set_title('Maximum Error')
    ax.set_xticks(range(len(tau_values)))
    ax.set_xticklabels([f'{v:.0e}' for v in tau_values], rotation=45)
    ax.set_yticks(range(len(discretization_factors)))
    ax.set_yticklabels([f'{v:.1f}' for v in discretization_factors])
    plt.colorbar(im1, ax=ax, label='Max Error')
    
    # Plot 2: Solve time heatmap
    ax = axes[0, 1]
    im2 = ax.imshow(solve_times, aspect='auto', cmap='plasma', origin='lower')
    ax.set_xlabel('tau index')
    ax.set_ylabel('discretization index')
    ax.set_title('Solve Time (seconds)')
    ax.set_xticks(range(len(tau_values)))
    ax.set_xticklabels([f'{v:.0e}' for v in tau_values], rotation=45)
    ax.set_yticks(range(len(discretization_factors)))
    ax.set_yticklabels([f'{v:.1f}' for v in discretization_factors])
    plt.colorbar(im2, ax=ax, label='Time (s)')

    # Plot 3: Steps heatmap (if available)
    if has_steps:
        ax = axes[0, 2]
        im3 = ax.imshow(num_steps, aspect='auto', cmap='cividis', origin='lower')
        ax.set_xlabel('tau index')
        ax.set_ylabel('discretization index')
        ax.set_title('Number of Internal Steps')
        ax.set_xticks(range(len(tau_values)))
        ax.set_xticklabels([f'{v:.0e}' for v in tau_values], rotation=45)
        ax.set_yticks(range(len(discretization_factors)))
        ax.set_yticklabels([f'{v:.1f}' for v in discretization_factors])
        plt.colorbar(im3, ax=ax, label='# steps')
    
    # Plot 3: Max error vs parameters (log scale)
    ax = axes[1, 0]
    for i, disc in enumerate(discretization_factors):
        ax.semilogy(tau_values, max_errors[i, :], 'o-', label=f'disc={disc:.1f}')
    ax.set_xlabel('tau_reaction')
    ax.set_ylabel('Max Error (log scale)')
    ax.set_title('Max Error vs tau (by discretization)')
    ax.set_xscale('log')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    
    # Plot 4: Solve time vs parameters
    ax = axes[1, 1]
    col_idx = 1 if not has_steps else 1
    ax = axes[1, col_idx]
    for i, disc in enumerate(discretization_factors):
        ax.plot(tau_values, solve_times[i, :], 'o-', label=f'disc={disc:.1f}')
    ax.set_xlabel('tau_reaction')
    ax.set_ylabel('Solve Time (seconds)')
    ax.set_title('Solve Time vs tau (by discretization)')
    ax.set_xscale('log')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Plot 6: Steps vs tau (if available)
    if has_steps:
        ax = axes[1, 2]
        for i, disc in enumerate(discretization_factors):
            ax.plot(tau_values, num_steps[i, :], 'o-', label=f'disc={disc:.1f}')
        ax.set_xlabel('tau_reaction')
        ax.set_ylabel('Number of Internal Steps')
        ax.set_title('Steps vs tau (by discretization)')
        ax.set_xscale('log')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save the figure
    # Save to configured results directory as PDF
    os.makedirs(path_to_save, exist_ok=True)
    output_file = os.path.join(path_to_save, os.path.basename(npz_file).replace('.npz', '_plots.pdf'))
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved plots to: {output_file}")

    # Only show plots if configured
    if show_plots:
        plt.show()
    else:
        plt.close('all')
    
    # Print summary statistics
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    
    # Count successful runs
    n_successful = np.sum(~np.isnan(max_errors))
    n_total = max_errors.size
    
    print(f"\nNumber of successful runs: {n_successful}/{n_total}")
    print(f"Number of failed/timed-out runs: {n_total - n_successful}/{n_total}")
    
    if n_successful == 0:
        print("\n⚠️  WARNING: All runs failed or timed out!")
        print("No valid results to analyze.")
        return
    
    # Find best parameters (minimum error)
    try:
        min_error_idx = np.unravel_index(np.nanargmin(max_errors), max_errors.shape)
        min_error = max_errors[min_error_idx]
        best_disc = discretization_factors[min_error_idx[0]]
        best_tau = tau_values[min_error_idx[1]]
        
        print(f"\nBest parameters (minimum error):")
        print(f"  discretization_factor = {best_disc:.2f}")
        print(f"  tau_reaction = {best_tau:.2e}")
        print(f"  Max error = {min_error:.6e}")
        print(f"  Solve time = {solve_times[min_error_idx]:.2f} s")
    except ValueError:
        print("\nBest parameters (minimum error): N/A (all runs failed)")
    
    # Find fastest parameters
    try:
        min_time_idx = np.unravel_index(np.nanargmin(solve_times), solve_times.shape)
        min_time = solve_times[min_time_idx]
        fastest_disc = discretization_factors[min_time_idx[0]]
        fastest_tau = tau_values[min_time_idx[1]]
        
        print(f"\nFastest parameters:")
        print(f"  discretization_factor = {fastest_disc:.2f}")
        print(f"  tau_reaction = {fastest_tau:.2e}")
        print(f"  Max error = {max_errors[min_time_idx]:.6e}")
        print(f"  Solve time = {min_time:.2f} s")
    except ValueError:
        print("\nFastest parameters: N/A (all runs failed)")
    
    # Overall statistics
    print(f"\nOverall statistics:")
    if n_successful > 0:
        print(f"  Max error range: [{np.nanmin(max_errors):.6e}, {np.nanmax(max_errors):.6e}]")
        print(f"  Solve time range: [{np.nanmin(solve_times):.2f} s, {np.nanmax(solve_times):.2f} s]")
        if has_steps:
            print(f"  Steps range: [{np.nanmin(num_steps):.0f}, {np.nanmax(num_steps):.0f}]")
    else:
        print(f"  Max error range: N/A")
        print(f"  Solve time range: N/A")


if __name__ == "__main__":
    load_and_plot_results("data/hyperparameter_search_ADAMS.npz")

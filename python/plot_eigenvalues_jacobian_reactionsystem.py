"""Plot eigenvalues of the reaction system Jacobian and a heatmap of deviation
from a heuristic target eigenvalue over (k_f_W, k_f_Tris).

First figure: real parts of numerical vs analytical eigenvalues.
Second figure: for each concentration set, heatmap of max(|lambda_i - lambda_target|)
for the two non-zero eigenvalues.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import linalg
from matplotlib.colors import LogNorm, Normalize, PowerNorm
from shared_data import path_to_save, show_plots

def analytical_lambda12(c, k_f_W, k_f_Tris, K_W, K_Tris):
    """Compute analytical eigenvalues (lambda1, lambda2) for reduced 2x2 system.
    Returns complex values to handle negative discriminants.
    """
    k_r_W = k_f_W / K_W
    k_r_Tris = k_f_Tris / K_Tris
    c_H2O, c_H_plus, c_OH_minus, c_TrisH_plus, c_Tris = c
    MN = np.array([
        [-(k_f_W + k_r_W*(c_H_plus + c_OH_minus)),   -k_r_W*c_OH_minus],
        [ (k_f_Tris - k_r_Tris*(c_H_plus + c_Tris)),  -(k_r_Tris*c_Tris + k_f_Tris)]
    ])
    tr = np.trace(MN)
    det = np.linalg.det(MN)
    discriminant = tr**2 - 4*det
    sqrt_discriminant = np.sqrt(discriminant + 0j)
    lambda1 = 0.5 * (tr + sqrt_discriminant)
    lambda2 = 0.5 * (tr - sqrt_discriminant)
    return lambda1, lambda2

def jacobian(c, k_f_W, k_r_W, k_f_Tris, k_b_Tris):
    """Return 5x5 Jacobian matrix for reactions:
    H2O <-> H+ + OH- and TrisH+ <-> H+ + Tris.
    c = (c_H2O, c_H_plus, c_OH_minus, c_TrisH_plus, c_Tris)
    """
    c_H2O, c_H_plus, c_OH_minus, c_TrisH_plus, c_Tris = c
    return np.array([
        [-k_f_W,  k_r_W*c_OH_minus,                        k_r_W*c_H_plus,              0,                     0],
        [ k_f_W, -k_r_W*c_OH_minus - k_b_Tris*c_Tris,      -k_r_W*c_H_plus,             k_f_Tris,             -k_b_Tris*c_H_plus],
        [ k_f_W, -k_r_W*c_OH_minus,                       -k_r_W*c_H_plus,              0,                     0],
        [ 0,      k_b_Tris*c_Tris,                         0,                           -k_f_Tris,             k_b_Tris*c_H_plus],
        [ 0,     -k_b_Tris*c_Tris,                         0,                            k_f_Tris,            -k_b_Tris*c_H_plus]
    ])

# Example equilibrium concentrations (arbitrary positive numbers)
concentration_arrays = [
    ("Initial concentrations", np.array([1.0, 1e-7, 1e-7, 0.0, 0.1])),
    ("Equilibrium concentrations", np.array([1.0, 10**-10.4, 10**-3.6, 0.1 * 10**-2.3, 0.1])),
    ("Other off-equilibrium concentrations", np.array([1.0, 0, 0.1, 0.1, 0.0]))
]

K_W = 10**-14
K_Tris = 10**-8.1

# Heuristic target eigenvalue (can be adjusted)
lambda_target = -1e2 # choose a representative negative rate; adjust as needed

# Global plotting style tweaks (slightly larger fonts and consistent math rendering)
plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 14,
    "legend.fontsize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "mathtext.default": "regular",
})

# Check rank and eigenvalues for all concentration arrays
for name, c in concentration_arrays:
    c_H2O, c_H_plus, c_OH_minus, c_TrisH_plus, c_Tris = c
    
    Q_W = (c_H_plus * c_OH_minus / c_H2O) if c_H2O != 0 else 0
    O_Tris = (c_H_plus * c_Tris / c_TrisH_plus) if c_TrisH_plus != 0 else 0
    
    print(f"\n{'='*60}")
    print(f"{name}:")
    print(f"pK_W = {-np.log10(K_W):.2f} vs pQ_W = {-np.log10(Q_W) if Q_W > 0 else 'inf'}")
    print(f"pK_Tris = {-np.log10(K_Tris):.2f} vs pQ_Tris = {-np.log10(O_Tris) if O_Tris > 0 else 'inf'}")

# Grids for first figure (keep original behavior)
k_f_W_values_fig = np.logspace(-15, 0, 6)
k_f_Tris_values_fig = np.logspace(-10, 5, 200)

# Grids for heatmap (denser in k_f_W, slightly less dense in k_f_Tris)
k_f_W_values_heat = np.logspace(-15, 0, 50)
k_f_Tris_values_heat = np.logspace(-10, 5, 50)

# Create subplots: rows for k_f_W values, columns for concentration arrays (first figure)
fig, axes = plt.subplots(len(k_f_W_values_fig), len(concentration_arrays),
                         figsize=(6*len(concentration_arrays), 4*len(k_f_W_values_fig)))

# Prepare containers for heatmap data for second figure
heatmaps = {}

# Ensure axes is 2D
if len(k_f_W_values_fig) == 1:
    axes = axes.reshape(1, -1)
if len(concentration_arrays) == 1:
    axes = axes.reshape(-1, 1)

for col_idx, (conc_name, c) in enumerate(concentration_arrays):
    for row_idx, k_f_W in enumerate(k_f_W_values_fig):
        k_r_W = k_f_W / K_W  # reverse rate constant for water autoionization
        
        numerical_eigen_real_parts = []
        analytical_eigen_real_parts = []
        
        for t_idx, k_f_Tris in enumerate(k_f_Tris_values_fig):
            k_r_Tris = k_f_Tris / K_Tris  # reverse rate constant for Tris dissociation

            J = jacobian(c, k_f_W, k_r_W, k_f_Tris, k_r_Tris)
            # Use scipy's eig with balancing for better numerical stability
            numerical_eigvals = linalg.eig(J, left=False, right=False)

            # Analytical eigenvalues (two non-zero of the full 5x5 Jacobian)
            lambda1, lambda2 = analytical_lambda12(c, k_f_W, k_f_Tris, K_W, K_Tris)

            analytical_eigvals = np.array([0, 0, 0, lambda1, lambda2])
            
            # Sort eigenvalues by real part
            numerical_eigvals_sorted = np.sort(numerical_eigvals.real)
            numerical_eigen_real_parts.append(numerical_eigvals_sorted)

            analytical_eigvals_sorted = np.sort(analytical_eigvals.real)
            analytical_eigen_real_parts.append(analytical_eigvals_sorted)
        
        numerical_eigen_real_parts = np.array(numerical_eigen_real_parts)
        analytical_eigen_real_parts = np.array(analytical_eigen_real_parts)
        
        # Plot eigenvalues (real parts)
        ax = axes[row_idx, col_idx]
        for i in range(numerical_eigen_real_parts.shape[1]):
            ax.plot(k_f_Tris_values_fig, numerical_eigen_real_parts[:, i], color='grey', alpha=0.7, linewidth=2)
            ax.plot(k_f_Tris_values_fig, analytical_eigen_real_parts[:, i], color='black', alpha=0.7, linestyle='--', linewidth=2)
        
        ax.set_xscale('log')
        ax.set_yscale('symlog', linthresh=1)
        ax.set_ylim(bottom=None, top=1)  # Only show negative values (y <= 0)
        
        # Only add x-label to bottom row
        if row_idx == len(k_f_W_values_fig) - 1:
            ax.set_xlabel(r"$k_{f,Tris}$")

        ax.set_ylabel(r"$\operatorname{Re}(\lambda)$")

        # Add title with concentration type and k_f_W value
        if row_idx == 0:
            ax.set_title(f"{conc_name}\n$k_{{f,W}}$ = {k_f_W:.2e}", pad=10)
        else:
            ax.set_title(f"$k_{{f,W}}$ = {k_f_W:.2e}", pad=10)

        ax.grid(True, which="both", ls="--", alpha=0.5)

    # First figure only

# Add a single legend for the entire figure
handles = [
    plt.Line2D([0], [0], color='grey', alpha=0.7, linewidth=2, label='Numerical'),
    plt.Line2D([0], [0], color='black', alpha=0.7, linewidth=2, linestyle='--', label='Analytical')
]
fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=True)

plt.tight_layout()
plt.subplots_adjust(top=0.95, hspace=0.3)  # Make room for legend and increase vertical spacing
os.makedirs(path_to_save, exist_ok=True)
plt.savefig(os.path.join(path_to_save, "jacobian_eigenvalues_H2O_Tris.pdf"), bbox_inches="tight")

# Build heatmap data on separate grids
eps = 1e-24
for conc_idx, (conc_name, c) in enumerate(concentration_arrays):
    heatmap = np.zeros((len(k_f_W_values_heat), len(k_f_Tris_values_heat)))
    # Stabilize lambda_target for log (absolute magnitude)
    for i_w, kfW in enumerate(k_f_W_values_heat):
        for i_t, kfT in enumerate(k_f_Tris_values_heat):
            l1, l2 = analytical_lambda12(c, kfW, kfT, K_W, K_Tris)
            # Distance in log-space of magnitudes: |log10(|Re(lambda)|) - log10(|lambda_target|)|
            diff1 = abs(np.log10(abs(l1.real) / abs(lambda_target) + 1e-20))
            diff2 = abs(np.log10(abs(l2.real) / abs(lambda_target) + 1e-20))
            heatmap[i_w, i_t] = max(diff1, diff2)
    heatmaps[conc_name] = heatmap

# Second figure: heatmaps of max(|lambda_i - lambda_target|)
fig_heat, axes2 = plt.subplots(1, len(concentration_arrays), figsize=(6*len(concentration_arrays), 5))
if len(concentration_arrays) == 1:
    axes2 = [axes2]

# Mesh for pcolormesh (centers) using heatmap grids
X, Y = np.meshgrid(k_f_Tris_values_heat, k_f_W_values_heat)
for conc_idx, (conc_name, _) in enumerate(concentration_arrays):
    data = heatmaps[conc_name]
    ax2 = axes2[conc_idx]
    # Emphasize low values smoothly using a power-law normalization (gamma < 1)
    vmax = float(np.nanmax(data)) if np.isfinite(np.nanmax(data)) else 1.0
    norm = PowerNorm(gamma=0.5, vmin=0.0, vmax=vmax)
    pcm = ax2.pcolormesh(X, Y, data, shading='auto', norm=norm, cmap='viridis')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    # Ensure full ranges are visible on both axes
    ax2.set_xlim(k_f_Tris_values_heat[0], k_f_Tris_values_heat[-1])
    ax2.set_ylim(k_f_W_values_heat[0], k_f_W_values_heat[-1])
    ax2.set_xlabel(r"$k_{f,Tris}$")
    ax2.set_ylabel(r"$k_{f,W}$")
    ax2.set_title(f"{conc_name}")
    cbar = fig_heat.colorbar(pcm, ax=ax2)
    cbar.set_label(r"$|\log_{10}(|\mathrm{Re}(\lambda_i)| / |\lambda_{\mathrm{target}}|)|$")

fig_heat.suptitle(rf"$|\log_{{10}}(\mathrm{{Re}}(\lambda_i) / \lambda_{{\mathrm{{target}}}} + 10^{{-20}})| \quad (\lambda_{{\mathrm{{target}}}}={lambda_target})$", y=0.98)
fig_heat.tight_layout()
fig_heat.savefig(os.path.join(path_to_save, "jacobian_eigenvalues_diff_heatmap.pdf"), bbox_inches="tight")

# Only show plots if configured
if show_plots:
    plt.show()
else:
    plt.close('all')

# Verification printout for specified rates across concentration sets
test_kfW = 1e-12
test_kfT = 10.0
print("\nVerification at k_f_W=1e-12, k_f_Tris=10:")
for name, c in concentration_arrays:
    l1, l2 = analytical_lambda12(c, test_kfW, test_kfT, K_W, K_Tris)
    print(f"- {name}:")
    print(f"  Analytical: lambda1={l1}, lambda2={l2}")

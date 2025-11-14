#!/usr/bin/env python3
"""
Plot f(σ0; n0) for several electrolyte concentrations with fixed defaults.

- No argument parser
- No file saving
- Shows the plot interactively
- All plotted lines use gray/black tones (no colors)

Defaults match the original script's no-argument run:
    ε_r = 78.5 (water), T = 298.15 K,
    σ0 ∈ [0.0, 0.2] with 400 steps,
    concentrations = [1e-4, 1e-3, 1e-2, 1e-1] mol/L.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
import os
from shared_data import path_to_save, show_plots

try:
    from scipy.constants import epsilon_0, k as k_B, Avogadro, e
except Exception as exc:  # pragma: no cover
    raise SystemExit(
        "This script requires scipy. Please install scipy (and numpy, matplotlib)."
    ) from exc


def molar_to_number_density(c_mol_per_L: float | np.ndarray) -> np.ndarray:
    """Convert concentration from mol/L to number density [1/m^3].

    1 L = 1e-3 m^3, thus 1 mol/L = 1000 mol/m^3.
    n0 = c [mol/m^3] * N_A.
    """
    c_mol_per_m3 = np.asarray(c_mol_per_L) * 1000.0
    return c_mol_per_m3 * Avogadro


def f_from_sigma(
    sigma0: np.ndarray,
    n0_number_density: float,
    epsilon_r: float = 78.5,
    T: float = 298.15,
) -> np.ndarray:
    """Compute f for given surface charge density array and number density.

    Parameters
    - sigma0: surface charge density [C/m^2]
    - n0_number_density: number density [1/m^3]
    - epsilon_r: relative permittivity of the medium [-]
    - T: temperature [K]
    """
    if n0_number_density <= 0:
        raise ValueError("n0 must be positive")

    eps = epsilon_r * epsilon_0
    denom = np.sqrt(2.0 * eps * k_B * T * n0_number_density)

    # From the provided equation
    a = -0.5 * sigma0 / denom
    b = np.sqrt((sigma0 ** 2) / (8.0 * eps * k_B * T * n0_number_density) + 1.0)
    f_alternative = (a + b) ** 2

    normalized_sigma0 = 0.5 * sigma0 / denom
    # Alternative expression for f, should be equivalent
    f = (-normalized_sigma0 + np.sqrt(normalized_sigma0 ** 2 + 1.0)) ** 2

    assert np.allclose(f_alternative, f), "f expressions do not match!"
    return f


if __name__ == "__main__":
    # Global plotting style tweaks for readability
    plt.rcParams.update({
        "font.size": 12,
        "axes.titlesize": 14,
        "axes.labelsize": 14,
        "legend.fontsize": 12,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "mathtext.default": "regular",
    })
    MNP_Ns = 1.08 * 1e-5
    maximum_possible_sigma = 1 * MNP_Ns * e * Avogadro

    sigma_min = -maximum_possible_sigma
    sigma_max = maximum_possible_sigma
    sigma_steps = 400
    epsilon_r = 78.5
    T = 298.15
    concentrations_mol_per_L = [1e-4, 1e-3, 1e-2, 1e-1, 1]

    sigma0 = np.linspace(sigma_min, sigma_max, sigma_steps)

    # Use seaborn-style grid but enforce grayscale for plotted curves
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(12, 8), constrained_layout=True)

    # Define grayscale tones for the four series (dark to light)
    gray_colors = ["#BABABA", "#828282", "#5D5D5D", "#3b3b3b","#000000"]

    for i, c in enumerate(concentrations_mol_per_L):
        n0 = float(molar_to_number_density(c))
        f_vals = f_from_sigma(sigma0, n0, epsilon_r=epsilon_r, T=T)
        # Legend in LaTeX with escaped braces; avoid Unicode superscripts
        label = rf"$c={c:g}\,\mathrm{{mol\,L^{{-1}}}}$ ($n_0={n0:.2e}\,\mathrm{{m^{{-3}}}}$)"
        color = gray_colors[i % len(gray_colors)]
        ax.plot(sigma0, f_vals, label=label, color=color, lw=2.0)

    ax.set_xlabel(r"$\sigma_0$ [C m$^{-2}$]")
    ax.set_ylabel(r"$f$")
    # ax.set_ylim(0,5)
    ax.set_yscale("log")
    ax.set_title(fr"$\varepsilon_r={epsilon_r}$, $T={T}\,K$")
    ax.legend(title="Electrolyte concentration")
    ax.axhline(1.0, color="k", lw=0.8, ls=":", alpha=0.7)
    ax.grid(True, which="both", alpha=0.25)

    # Save to configured results directory as PDF
    os.makedirs(path_to_save, exist_ok=True)
    outfile = os.path.join(path_to_save, "gouy_chapman_f_vs_sigma.pdf")
    plt.savefig(outfile, bbox_inches="tight")

    # Only show plots if configured
    if show_plots:
        plt.show()
    else:
        plt.close('all')

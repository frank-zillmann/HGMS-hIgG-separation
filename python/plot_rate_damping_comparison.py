#!/usr/bin/env python3
"""Compare reaction rate damping strategies.

Generates a plot of the physical reaction rate (identity) vs three damping
functions:
  1. Hard clip (Eq. \ref{Eq:rate_clip})
  2. Smooth (Eq. \ref{Eq:rate_smooth})
  3. Log compression (Eq. \ref{Eq:rate_log})

Adjust nu_max or alpha below if desired.

Run:
    python plot_rate_damping_comparison.py

Creates: rate_damping_comparison.png in current directory.
"""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
import os
from shared_data import path_to_save, show_plots
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

nu_max = 1.0       # heuristic maximum rate magnitude
alpha = 1.0        # scaling for log compression
n_points = 1000

nu_phys = np.linspace(-5*nu_max, 5*nu_max, n_points)

# Hard clip
nu_clip = np.clip(nu_phys, -nu_max, nu_max)

# Smooth rational saturation (Eq. rate_smooth)
nu_smooth = nu_phys / (1.0 + np.abs(nu_phys)/nu_max)

# Log compression (Eq. rate_log)
nu_log = np.sign(nu_phys) * (1.0/alpha) * np.log(1.0 + alpha * np.abs(nu_phys))

fig, ax = plt.subplots(figsize=(12,8))
ax.plot(nu_phys, nu_phys, label="physical (identity)", color="0.0", linewidth=2.0, linestyle="-")
ax.plot(nu_phys, nu_clip, label="hard clipping", color="0.3", linewidth=2.0, linestyle="--")
ax.plot(nu_phys, nu_smooth, label="smooth limiting", color="0.5", linewidth=2.0, linestyle="-.")
ax.plot(nu_phys, nu_log, label="log compression", color="0.7", linewidth=2.0, linestyle=":")

ax.axhline(nu_max, color="0.7", linestyle="--", linewidth=1.0)
ax.axhline(-nu_max, color="0.7", linestyle="--", linewidth=1.0)
ax.text(0.02, nu_max+0.03, r"$+\nu_{max}$", color="0.3")
ax.text(0.02, -nu_max-0.09, r"$-\nu_{max}$", color="0.3")

ax.set_xlabel(r"physical rate $\nu$")
ax.set_ylabel(r"damped rate $\tilde{\nu}$")
ax.set_title(r"Reaction rate damping strategies")
ax.legend(frameon=False)
ax.grid(alpha=0.3)

fig.tight_layout()

os.makedirs(path_to_save, exist_ok=True)
fig.savefig(os.path.join(path_to_save, "rate_damping_comparison.pdf"), bbox_inches="tight")
if __name__ == "__main__":
  if show_plots:
    plt.show()
  else:
    plt.close('all')

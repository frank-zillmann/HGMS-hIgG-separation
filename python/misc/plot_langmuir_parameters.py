import matplotlib.pyplot as plt
import numpy as np
import os
from shared_data import path_to_save, show_plots


# _fit data
pH_fit = np.array([7.0, 4.5, 3.7, 3.5, 2.9, 2.6, 2.5, 2.0])
K_B_fit = np.array([50478.7, 20176.3, 576.3, 484.8, 462.8, 342.8, 320, 5.2])/600
qMax_fit = np.array([0.3563, 0.2863, 0.134963, 0.05562, 0.05462, 0.024, 0.00001, 0])

# _fit_paper data
pH_fit_paper = np.array([7.0, 4.5, 3.7, 3.5, 2.6, 2.5, 2.0])
K_B_fit_paper = np.array([36.52, 35.15, 3.84, 3.22, 2.28, 2.14, 0.03])
qMax_fit_paper = np.array([0.31, 0.28, 0.13, 0.05, 0.02, 0.01, 0.0])

# _exp data
pH_exp = np.array([7.5, 7.0, 4.5, 3.5, 2.9])
K_B_exp = np.array([42.63, 95.98, 93.5, 53.52, 53.28])
qMax_exp = np.array([0.15, 0.13, 0.18, 0.03, 0.04])

plt.figure(figsize=(10, 10))

plt.subplot(1, 2, 1)
plt.plot(pH_fit, K_B_fit, marker='o', color='orange', label='K_B_fit')
plt.plot(pH_fit_paper, K_B_fit_paper, marker='s', color='green', label='K_B_fit_paper')
plt.plot(pH_exp, K_B_exp, marker='^', color='red', label='K_B_exp')
plt.xlabel('pH')
plt.ylabel('K_B')
plt.title('K_B vs pH')
plt.legend()

from scipy.optimize import curve_fit

def atan_sigmoid(x, q_inf, pH_crit, pH_width):
	return q_inf/2 * ((2/np.pi) * np.arctan(pH_width*(x - pH_crit)) + 1)

 # Fit the curve to the _fit data
popt_fit, _ = curve_fit(atan_sigmoid, pH_fit, qMax_fit, p0=[0.36, 3.5, 1])
q_inf_fit, pH_crit_fit, pH_width_fit = popt_fit

plt.subplot(1, 2, 2)
plt.plot(pH_fit, qMax_fit, marker='o', color='orange', label='qMax_fit')
plt.plot(pH_fit_paper, qMax_fit_paper, marker='s', color='green', label='qMax_fit_paper')
plt.plot(pH_exp, qMax_exp, marker='^', color='red', label='qMax_exp')
# Plot fitted curve for _fit data
pH_curve = np.linspace(min(pH_fit), max(pH_fit), 200)
plt.plot(pH_curve, atan_sigmoid(pH_curve, q_inf_fit, pH_crit_fit, pH_width_fit), color='blue', label=f'Fit: q_inf={q_inf_fit:.3f}, pH_crit={pH_crit_fit:.3f}, pH_width={pH_width_fit:.3f}')
plt.xlabel('pH')
plt.ylabel('qMax')
plt.title('qMax vs pH')
plt.legend()

plt.tight_layout()

# Save to configured results directory as PDF
os.makedirs(path_to_save, exist_ok=True)
outfile = os.path.join(path_to_save, "plot_langmuir_parameters.pdf")
plt.savefig(outfile, bbox_inches="tight")

# Only show plots if configured
if show_plots:
	plt.show()
else:
	plt.close('all')
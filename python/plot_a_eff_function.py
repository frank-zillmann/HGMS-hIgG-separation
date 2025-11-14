import numpy as np
import matplotlib.pyplot as plt
import os
from shared_data import path_to_save, show_plots

def a_eff_function(um_u0_ratio):
    a1 = 2.035
    a2 = 107.1
    a3 = -0.00808
    p = 0.07477
    q = 1.083
    return 1 / (1 + a1 * np.exp(- p * um_u0_ratio) + a2 * np.exp(- q * um_u0_ratio) + a3)

x = np.linspace(0, 100, 500)
y = a_eff_function(x)

print(f"Left boundary: x={x[0]}, y={y[0]}")
print(f"Right boundary: x={x[-1]}, y={y[-1]}")

plt.figure(figsize=(8, 5))
plt.plot(x, y, label=r'$a_{eff}(um/u_0)$')
plt.xlabel(r'$um/u_0$')
plt.ylabel(r'$a_{eff}$')
plt.title('Effective Capture Area Function')
plt.grid(True)
plt.legend()
plt.tight_layout()

# Save to configured results directory as PDF
os.makedirs(path_to_save, exist_ok=True)
outfile = os.path.join(path_to_save, "plot_a_eff_function.pdf")
plt.savefig(outfile, bbox_inches="tight")

# Only show plots if configured
if show_plots:
    plt.show()
else:
    plt.close('all')

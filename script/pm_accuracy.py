import matplotlib.pyplot as plt
from load_data import *
import numpy as np
import os
import sys

OPTION_SOFT = 'soft'
OPTION_NO_SOFT = 'no-soft'

if len(sys.argv) != 2:
    print('Missing argument: soft / no-soft')
    sys.exit(1)

option = sys.argv[1]
if option not in [OPTION_SOFT, OPTION_NO_SOFT]:
    print('Unknown argument. Available are soft / no-soft')

# Configuration
schemes = ['NGP', 'CIC', 'TSC']
fd_orders = ['2-point', '4-point']  # 2-point and 4-point
grid_sizes = range(30, 70, 4)

# Plot setup
plt.figure()

if option == OPTION_NO_SOFT:
    forcesPP = load_data(rf"output\force-pp-0.dat")

# Iterate over interpolation schemes and finite difference orders
for scheme in schemes:
    for order in fd_orders:
        mean_deltas_pm = []
        valid_grid_sizes = []

        for i in grid_sizes:
            filename = rf"output\force-pm-{scheme}-{order}-{i}.dat"
            if not os.path.isfile(filename):
                print(filename, 'not found!')

            forcesPM = load_data(filename)
            if option == OPTION_SOFT:
                forcesPP = load_data(rf"output\force-pp-{i}.dat")
            deltas = []
            for f1, f2 in zip(forcesPP, forcesPM):
                num = np.sqrt((f1[0] - f2[0])**2 + (f1[1] - f2[1])**2 + (f1[2] - f2[2])**2)
                denom = np.sqrt(f1[0]**2 + f1[1]**2 + f1[2]**2)
                deltas.append(num / (denom + 1e-10))
            mean_deltas_pm.append(np.mean(deltas))
            valid_grid_sizes.append(i)

        label = f"{scheme.upper()}-{order}-point"
        plt.plot(valid_grid_sizes, mean_deltas_pm, marker=('o'if order == '4-point' else 's'), linestyle=('-' if order=='4-point' else '--'), label=label)

# Final plot adjustments
plt.xlabel(r'Grid size')
plt.ylabel('Mean Relative Error')
plt.title('PP with' + ('' if option == OPTION_SOFT else 'out') + ' Force Softening')
plt.yscale('log')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

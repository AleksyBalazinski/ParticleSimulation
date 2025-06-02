import matplotlib.pyplot as plt
from load_data import *
import numpy as np
import os

# Reference force data (particle-particle method)
forcesPP = load_data(r"output\force-pp.dat")

# Configuration
schemes = ['TSC', 'CIC', 'NGP']
fd_orders = ['2-point', '4-point']  # 2-point and 4-point
grid_sizes = range(5, 25, 1)

# Plot setup
#plt.figure(figsize=(10, 6))

# Iterate over interpolation schemes and finite difference orders
for scheme in schemes:
    for order in fd_orders:
        mean_deltas_pm = []
        valid_grid_sizes = []

        for i in grid_sizes:
            filename = rf"output\force-p3m-{scheme}-{order}-{i}.dat"
            #filename = rf"build\source\Release\force-pm-lap-tsc-4point-64.dat"
            if not os.path.isfile(filename):
                print(filename, 'not found')
                continue  # Skip missing files

            forcesPM = load_data(filename)
            deltas = []
            for f1, f2 in zip(forcesPP, forcesPM):
                num = np.sqrt((f1[0] - f2[0])**2 + (f1[1] - f2[1])**2 + (f1[2] - f2[2])**2)
                denom = np.sqrt(f1[0]**2 + f1[1]**2 + f1[2]**2)
                deltas.append(num / (denom + 1e-10))
            mean_deltas_pm.append(np.mean(deltas))
            valid_grid_sizes.append(0.2 * i)

        label = f"{scheme.upper()}-{order}-point"
        plt.plot(valid_grid_sizes, mean_deltas_pm, marker=('o' if order=='4-point' else 's'), linestyle=('-' if order=='4-point' else '--'), label=label)

# Final plot adjustments
plt.xlabel(r'$a/H$')
plt.ylabel('Mean Relative Error')
plt.yscale('log')
plt.title('Error vs Particle Diameter')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

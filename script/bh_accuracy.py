import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import sys
from load_data import *
import numpy as np

forcesPP = load_data(r"output\force-pp.dat")
mean_deltas_quad = []
mean_deltas_mon = []
thetas = []

for i in range(12):
    thetas.append(i * 0.1)
    forcesBH_mon = load_data(rf"output\force-bh-{i}-m.dat")
    forcesBH_quad = load_data(rf"output\force-bh-{i}-q.dat")
    deltas = []
    for f1, f2 in zip(forcesPP, forcesBH_mon):
        num = np.sqrt((f1[0] - f2[0])**2 + (f1[1] - f2[1])**2 + (f1[2] - f2[2])**2)
        denom = np.sqrt(f1[0]**2 + f1[1]**2 + f1[2]**2)
        deltas.append(num / (denom + 1e-10))
    mean_deltas_mon.append(np.mean(deltas))

    forcesBH_mon = load_data(rf"output\force-bh-{i}-m.dat")
    deltas = []
    for f1, f2 in zip(forcesPP, forcesBH_quad):
        num = np.sqrt((f1[0] - f2[0])**2 + (f1[1] - f2[1])**2 + (f1[2] - f2[2])**2)
        denom = np.sqrt(f1[0]**2 + f1[1]**2 + f1[2]**2)
        deltas.append(num / (denom + 1e-10))
    mean_deltas_quad.append(np.mean(deltas))

plt.plot(thetas, mean_deltas_mon, marker='s', linestyle='-', color='r', label='Monopole Approximation')
plt.plot(thetas, mean_deltas_quad, marker='o', linestyle='-', color='b', label='Quadrupole Approximation')
plt.xlabel(r'$\theta$')
plt.ylabel('Mean Relative Error')
plt.yscale('log')
plt.title('Error vs Opening Angle')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
import matplotlib.pyplot as plt
import numpy as np
from load_data import *

m = 1
G = 4.5e-3

def s1_smoothed_g(rs, a):
    gs = []
    for r in rs:
        if r <= a:
            gs.append(G * m / a**2 * (8 * r / a - 9 * r**2 / a**2 + 2 * r**4 / a**4))
        else:
            gs.append(G * m / r**2)
    return gs

def s2_smoothed_g(rs, a):
    gs = []
    for r in rs:
        u = 2 * r / a
        if u <= 1:
            gs.append(G * m / (35 * a**2) * (224 * u - 224 * u**3 + 70 * u**4 + 48 * u**5 - 21 * u**6))
        elif u <= 2:
            gs.append(G * m / (35 * a**2) * (12 / u**2 - 224 + 896 * u - 840 * u**2 + 224 * u**3 + 70 * u**4 - 48 * u**5 + 7 * u**6))
        else:
            gs.append(G * m / r**2)
    return gs

def mean_normalized_rmse(y_true, y_pred):
    return np.sqrt(np.mean((y_true - y_pred) ** 2)) / np.mean(y_true) * 100

deviation_s1 = []
deviation_s2 = []
a_over_H_values = []

for i in range(5, 31):
    H = 1.875
    a = 0.2 * i * H
    a_over_H_values.append(0.2 * i)

    # S1
    frames_s1 = load_data(rf'output/positions-s1.dat')
    field_s1 = load_data(rf'output/field-s1-{i}.dat')

    x_s1 = frames_s1[-1][:, 0]
    g_x_s1 = field_s1[-1][:, 0]
    r_s1 = x_s1 - 30
    smoothed_g_x_s1 = s1_smoothed_g(r_s1, a)

    deviation_s1.append(mean_normalized_rmse(smoothed_g_x_s1, -np.array(g_x_s1)))

    # S2
    frames_s2 = load_data(rf'output/positions-s2.dat')
    field_s2 = load_data(rf'output/field-s2-{i}.dat')

    x_s2 = frames_s2[-1][:, 0]
    g_x_s2 = field_s2[-1][:, 0]
    r_s2 = x_s2 - 30
    smoothed_g_x_s2 = s2_smoothed_g(r_s2, a)

    deviation_s2.append(mean_normalized_rmse(smoothed_g_x_s2, -np.array(g_x_s2)))

    if i == 20:
        r_trunc = r_s1[5:]
        expected_g_x = G * m/r_trunc**2
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlim(0, 6)
        ax.set_ylim(0, 0.08)
        ax.set_title(f"$a = {int(0.2 * i)}H$")

        ax.set_xlabel('$r/H$')
        ax.set_ylabel('$F / (Gm_1m_2)$')


        ax.plot(np.array(r_s1)/H, -np.array(g_x_s1)/G, label='PM $S_1$', linestyle="-", color="red", linewidth=2)
        ax.plot(np.array(r_s1) / H, np.array(smoothed_g_x_s1)/G, label='Reference $S_1$', linestyle="--", color="c", linewidth=2)
        ax.plot(np.array(r_s2)/H, -np.array(g_x_s2)/G, label='PM $S_2$', linestyle="-", color="blue", linewidth=2)
        ax.plot(np.array(r_s2) / H, np.array(smoothed_g_x_s2)/G, label='Reference $S_2$', linestyle="--", color="green", linewidth=2)
        ax.plot(np.array(r_trunc) / H, np.array(expected_g_x)/G, label='Inverse square', linestyle=":", color="red", linewidth=2)
        ax.legend()

        plt.show()

# Plotting
plt.figure()
plt.plot(a_over_H_values, deviation_s1, marker='o', label='$S_1$')
plt.plot(a_over_H_values, deviation_s2, marker='s', label='$S_2$')
plt.xlabel("$a/H$")
plt.ylabel("Normalized RMSD (%)")
plt.yscale('log')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.show()

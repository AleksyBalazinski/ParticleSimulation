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

def process_frame_data(frames, fields, center):
    last_frame = frames[-1]
    last_field = fields[-1]

    x = last_frame
    g_x = last_field[:, 0]
    g_y = last_field[:, 1]
    g_z = last_field[:, 2]
    g = np.sqrt(g_x**2 + g_y**2 + g_z**2)

    r = np.sqrt((x[:, 0] - center[0])**2 + (x[:, 1] - center[1])**2 + (x[:, 2] - center[2])**2)
    return r, g

# Load data

frames_s1 = load_data(r'output/positions-s1.dat')
field_s1_ngp = load_data(r'output/field-s1-NGP.dat')
field_s1_cic = load_data(r'output/field-s1-CIC.dat')
field_s1_tsc = load_data(r'output/field-s1-TSC.dat')

# Define smoothing parameters
H = 1.875
a = 4 * H

# Process datasets
r_s1, g_s1_ngp = process_frame_data(frames_s1, field_s1_ngp, center=(30, 30, 15))
r_s1, g_s1_cic = process_frame_data(frames_s1, field_s1_cic, center=(30, 30, 15))
r_s1, g_s1_tsc = process_frame_data(frames_s1, field_s1_tsc, center=(30, 30, 15))
# Compute expected and smoothed values for original
expected_g_x = G * m / r_s1**2
s1_g_x = s1_smoothed_g(r_s1, a)

# Plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(0, 6)
ax.set_ylim(0, 0.1)

ax.set_xlabel('$r/H$')
ax.set_ylabel('$F / (Gm_1m_2)$')

ax.plot(np.array(r_s1)/H, np.array(g_s1_ngp)/G, label='NGP', linestyle="-", color="green", linewidth=1.5)
ax.plot(np.array(r_s1)/H, np.array(g_s1_cic)/G, label='CIC', linestyle="-", color="black", linewidth=1.5)
ax.plot(np.array(r_s1)/H, np.array(g_s1_tsc)/G, label='TSC', linestyle="-", color="black", linewidth=1.5)
ax.plot(np.array(r_s1[2:])/H, np.array(expected_g_x[2:])/G, label='Inverse square', linestyle=":", color="red", linewidth=1.5)
ax.plot(np.array(r_s1)/H, np.array(s1_g_x)/G, label='S1', linestyle="--", color="blue", linewidth=1.5)
ax.set_title(r'$a = 4H$')

ax.legend()
plt.show()
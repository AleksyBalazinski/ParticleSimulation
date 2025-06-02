import matplotlib.pyplot as plt
import numpy as np
from load_data import *

m = 1
G = 4.5e-3

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
frames_pm_tsc = load_data(r'output/positions-tsc.dat')
field_pm_tsc = load_data(r'output/field-tsc.dat')

frames_pm_cic = load_data(r'output/positions-cic.dat')
field_pm_cic = load_data(r'output/field-cic.dat')

frames_pm_ngp = load_data(r'output/positions-ngp.dat')
field_pm_ngp = load_data(r'output/field-ngp.dat')

# Define smoothing parameters
H = 1.875
a = 4 * H

# Process datasets
r_pm_tsc, g_pm_tsc = process_frame_data(frames_pm_tsc, field_pm_tsc, center=(30, 30, 15))
r_pm_cic, g_pm_cic = process_frame_data(frames_pm_cic, field_pm_cic, center=(30, 30, 15))
r_pm_ngp, g_pm_ngp = process_frame_data(frames_pm_cic, field_pm_ngp, center=(30, 30, 15))


# Compute expected and smoothed values for original
expected_g_x = G * m / r_pm_tsc**2

# Plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(0, 6)
ax.set_ylim(0, 0.5)

ax.set_xlabel('$r/H$')
ax.set_ylabel('$F / (Gm_1m_2)$')

ax.plot(np.array(r_pm_tsc)/H, np.array(g_pm_tsc)/G, label='TSC', linestyle="-", color="green", linewidth=1.5)
ax.plot(np.array(r_pm_cic)/H, np.array(g_pm_cic)/G, label='CIC', linestyle="-", color="black", linewidth=1.5)
ax.plot(np.array(r_pm_ngp)/H, np.array(g_pm_ngp)/G, label='NGP', linestyle="-", color="orange", linewidth=1.5)
ax.plot(np.array(r_pm_tsc[2:])/H, np.array(expected_g_x[2:])/G, label='Inverse square', linestyle=":", color="red", linewidth=1.5)

ax.legend()
plt.show()

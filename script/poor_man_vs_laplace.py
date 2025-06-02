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
frames_poor_man = load_data(r'output/positions-poor-man.dat')
field_poor_man = load_data(r'output/field-poor-man.dat')

frames_lap = load_data(r'output/positions-discrete-lap.dat')
field_lap = load_data(r'output/field-discrete-lap.dat')

# Define smoothing parameters
H = 1.875
a = 4 * H

# Process datasets
r_poor_man, g_poor_man = process_frame_data(frames_poor_man, field_poor_man, center=(30, 30, 15))
r_lap, g_lap = process_frame_data(frames_lap, field_lap, center=(30, 30, 15))

# Compute expected and smoothed values for original
expected_g_x = G * m / r_poor_man**2

# Plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(0, 6)
ax.set_ylim(0, 0.5)

ax.set_xlabel('$r/H$')
ax.set_ylabel('$F / (Gm_1m_2)$')

ax.plot(np.array(r_poor_man)/H, np.array(g_poor_man)/G, label="'poor man'", linestyle="-", color="green", linewidth=1.5)
ax.plot(np.array(r_lap)/H, np.array(g_lap)/G, label='discretized Laplacian', linestyle="-", color="black", linewidth=1.5)
ax.plot(np.array(r_lap[2:])/H, np.array(expected_g_x[2:])/G, label='Inverse square', linestyle=":", color="red", linewidth=1.5)

ax.legend()
plt.show()


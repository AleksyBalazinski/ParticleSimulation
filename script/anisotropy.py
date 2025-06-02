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
frames = load_data(r'output/positions.dat')
field = load_data(r'output/field.dat')

frames_offset = load_data(r'output/positions_offset.dat')
field_offset = load_data(r'output/field_offset.dat')

frames_diag = load_data(r'output/positions_diag.dat')
field_diag = load_data(r'output/field_diag.dat')

# Define smoothing parameters
H = 1.875
a = 4 * H

# Process datasets
r, g = process_frame_data(frames, field, center=(30, 30, 15))
r_offset, g_offset = process_frame_data(frames_offset, field_offset, center=(29, 30, 15))
r_diag, g_diag = process_frame_data(frames_diag, field_diag, center=(30, 30, 15))

# Compute expected and smoothed values for original
expected_g_x = G * m / r**2

# Plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(0, 10)
ax.set_ylim(0, 0.15)

ax.set_xlabel('$r/H$')
ax.set_ylabel('$F / (Gm_1m_2)$')

ax.plot(np.array(r)/H, np.array(g)/G, label='PM (along x-axis)', linestyle="-", color="blue", linewidth=1.5)
ax.plot(np.array(r_offset)/H, np.array(g_offset)/G, label='PM (offset source)', linestyle="-", color="green", linewidth=1.5)
ax.plot(np.array(r_diag)/H, np.array(g_diag)/G, label='PM (along $x=y$ diagonal)', linestyle="-", color="orange", linewidth=1.5)

ax.plot(np.array(r[5:])/H, np.array(expected_g_x[5:])/G, label='Inverse square', linestyle=":", color="red", linewidth=1.5)

ax.axvline(x=1, color='black', linestyle='--', linewidth=1.5, label='$r=H$')
ax.axvline(x=0.5, color='gray', linestyle='--', linewidth=1.5, label='$r=H/2$')

ax.legend()
plt.show()

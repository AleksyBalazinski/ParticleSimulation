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
frames_2 = load_data(r'output/positions-2.dat')
field_2 = load_data(r'output/field-2.dat')

frames_4 = load_data(r'output/positions-4.dat')
field_4 = load_data(r'output/field-4.dat')

# Define smoothing parameters
H = 1.875
a = 4 * H

# Process datasets
r_2, g_2 = process_frame_data(frames_2, field_2, center=(30, 30, 15))
r_4, g_4 = process_frame_data(frames_4, field_4, center=(30, 30, 15))


# Compute expected and smoothed values for original
expected_g_x = G * m / r_2**2

# Plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(0, 6)
ax.set_ylim(0, 0.1)

ax.set_xlabel('$r/H$')
ax.set_ylabel('$F / (Gm_1m_2)$')

ax.plot(np.array(r_2)/H, np.array(g_2)/G, label='2-point', linestyle="-", color="green", linewidth=1.5)
ax.plot(np.array(r_4)/H, np.array(g_4)/G, label='4-point', linestyle="-", color="black", linewidth=1.5)
ax.plot(np.array(r_2[2:])/H, np.array(expected_g_x[2:])/G, label='Inverse square', linestyle=":", color="red", linewidth=1.5)
ax.set_title(r'$a = 4H$')

ax.legend()
plt.show()

error_field2 = np.abs((g_2[20:] - expected_g_x[20:]) / expected_g_x[20:])
error_field4 = np.abs((g_4[20:] - expected_g_x[20:]) / expected_g_x[20:])

r_err = r_2[20:] / H  # Normalized r for plotting

# Plot error
plt.plot(r_err, error_field2)
plt.plot(r_err, error_field4)

plt.xlabel('$r/H$')
plt.ylabel('Relative Error')
plt.xlim(3, 15)
plt.title('Relative Error vs. $r/H$')
plt.legend()
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.show()

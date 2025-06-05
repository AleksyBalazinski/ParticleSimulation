import matplotlib.pyplot as plt
import numpy as np

# Load data
z_order_data = np.loadtxt('bh_z_order.txt')
standard_data = np.loadtxt('bh_standard.txt')

# Scale particle numbers by 10,000
scale_factor = 10000
z_order_particles = z_order_data[:, 0] / scale_factor
z_order_time = z_order_data[:, 1]
standard_particles = standard_data[:, 0] / scale_factor
standard_time = standard_data[:, 1]

# Plotting
plt.plot(z_order_particles, z_order_time, 'o-', label='Z-order BH', color='blue')
plt.plot(standard_particles, standard_time, 's-', label='Standard BH', color='green')

# Labels and title
plt.title('Tree Construction Time vs Particle Number')
plt.xlabel(r'Number of Particles ($\times 10^4$)')
plt.ylabel('Time (seconds)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

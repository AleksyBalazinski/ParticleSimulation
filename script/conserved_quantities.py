import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# Helper function to apply scientific notation
def apply_sci_notation(ax):
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-3, 3))  # Adjust range if needed
    formatter.set_useOffset(False)
    ax.yaxis.set_major_formatter(formatter)

# energy
energy_data = np.loadtxt("output/energy.txt")
ep = energy_data[:, 0]
ek = energy_data[:, 1]
e = ep + ek
time = np.arange(len(ep))

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
ax1.plot(time, ep, label='EP', color='red')
ax1.plot(time, ek, label='EK', color='green')
ax1.plot(time, e, label='E', color='blue')
ax1.set_xlabel('Time')
ax1.set_ylabel('Energy')
ax1.legend()
apply_sci_notation(ax1)

# angular momentum
angular_momentum_data = np.loadtxt("output/angular_momentum.txt")
l_x = angular_momentum_data[:, 0]
l_y = angular_momentum_data[:, 1]
l_z = angular_momentum_data[:, 2]

fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111)
ax2.plot(time, l_x, label='$l_x$', color='blue')
ax2.plot(time, l_y, label='$l_y$', color='green')
ax2.plot(time, l_z, label='$l_z$', color='red')
ax2.set_xlabel('Time')
ax2.set_ylabel('Angular momentum components')
ax2.legend()
apply_sci_notation(ax2)

# momentum
momentum_data = np.loadtxt("output/momentum.txt")
p_x = momentum_data[:, 0]
p_y = momentum_data[:, 1]
p_z = momentum_data[:, 2]

expected_momentum_data = np.loadtxt("output/expected_momentum.txt")
pe_x = expected_momentum_data[:, 0]
pe_y = expected_momentum_data[:, 1]
pe_z = expected_momentum_data[:, 2]
time = np.arange(len(p_x))

fig3 = plt.figure(3)
ax3 = fig3.add_subplot(111)
ax3.plot(time, p_x, label='$p_x$', color='blue')
ax3.plot(time, p_y, label='$p_y$', color='green')
ax3.plot(time, p_z, label='$p_z$', color='red')

ax3.plot(time, pe_x, label='$p_x$ (expected)', color='blue', linestyle='--')
ax3.plot(time, pe_y, label='$p_y$ (expected)', color='green', linestyle='--')
ax3.plot(time, pe_z, label='$p_z$ (expected)', color='red', linestyle='--')

ax3.set_xlabel('Time')
ax3.set_ylabel('Momentum components')
ax3.legend()
apply_sci_notation(ax3)

plt.show()

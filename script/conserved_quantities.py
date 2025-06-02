import matplotlib.pyplot as plt
import numpy as np
import sys

# energy
energy_data = np.loadtxt("output/energy.txt")
ep = energy_data[:, 0]
ek = energy_data[:, 1]
e = ep + ek
time = np.arange(len(ep))

plt.figure(1)
plt.plot(time, ep, label='EP', color='red')
plt.plot(time, ek, label='EK', color='green')
plt.plot(time, e, label='E', color='blue')
plt.xlabel('Time')
plt.ylabel('Energy')
plt.legend()

# angular momentum
angular_momentum_data = np.loadtxt("output/angular_momentum.txt")
l_x = angular_momentum_data[:, 0]
l_y = angular_momentum_data[:, 1]
l_z = angular_momentum_data[:, 2]
plt.figure(2)
plt.plot(time, l_x, label='$l_x$', color='blue')
plt.plot(time, l_y, label='$l_y$', color='green')
plt.plot(time, l_z, label='$l_z$', color='red')
plt.xlabel('Time')
plt.ylabel('Angular momentum components')
plt.legend()

# momentum
momentum_data = np.loadtxt("output/momentum.txt")
p_x = momentum_data[:, 0]
p_y = momentum_data[:, 1]
p_z = momentum_data[:, 2]

# expected momentum
expected_momentum_data = np.loadtxt("output/expected_momentum.txt")
pe_x = expected_momentum_data[:, 0]
pe_y = expected_momentum_data[:, 1]
pe_z = expected_momentum_data[:, 2]
time = np.arange(len(p_x))

plt.figure(3)
plt.plot(time, p_x, label='$p_x$', color='blue')
plt.plot(time, p_y, label='$p_y$', color='green')
plt.plot(time, p_z, label='$p_z$', color='red')
if(len(sys.argv) >= 5):
    plt.plot(time, pe_x, label='$p_x$ (expected)', color='blue', linestyle='--')
    plt.plot(time, pe_y, label='$p_y$ (expected)', color='green', linestyle='--')
    plt.plot(time, pe_z, label='$p_z$ (expected)', color='red', linestyle='--')

plt.xlabel('Time')
plt.ylabel('Momentum components')
plt.legend()

plt.show()
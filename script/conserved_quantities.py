import matplotlib.pyplot as plt
import numpy as np
import sys

# energy
energy_data = np.loadtxt(sys.argv[1])
e = energy_data[:]
time = np.arange(len(e))

plt.figure(1)
plt.plot(time, e, label='E', color='blue')
plt.xlabel('Time')
plt.ylabel('Energy')
plt.legend()

# momentum
momentum_data = np.loadtxt(sys.argv[2])
p_x = momentum_data[:, 0]
p_y = momentum_data[:, 1]
p_z = momentum_data[:, 2]
time = np.arange(len(p_x))

plt.figure(2)
plt.plot(time, p_x, label='p_x', color='blue')
plt.plot(time, p_y, label='p_y', color='green')
plt.plot(time, p_z, label='p_z', color='red')

plt.xlabel('Time')
plt.ylabel('Momentum components')
plt.legend()

plt.show()
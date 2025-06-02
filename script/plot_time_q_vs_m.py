import matplotlib.pyplot as plt

# Helper function to read data from file
def read_times(filename):
    thetas = []
    times = []
    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split(',')
            thetas.append(0.1 * int(parts[0]))
            times.append(int(parts[1]))
    return thetas, times

# Read both datasets
mono_indices, mono_times = read_times('monopole_times.txt')
quad_indices, quad_times = read_times('quadrupole_times.txt')

# Plotting
plt.plot(mono_indices, mono_times, label='Monopole', marker='s', color='red')
plt.plot(quad_indices, quad_times, label='Quadrupole', marker='o', color='blue')

plt.xlabel(r'$\theta$')
plt.ylabel('Execution Time (ms)')
plt.yscale('log')
plt.title('Execution Time vs Opening Angle')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

import matplotlib.pyplot as plt

# Helper function to read data from file
def read_times(filename):
    thetas = []
    times = []
    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split(',')
            thetas.append(float(parts[0]))
            times.append(int(parts[1]))
    return thetas, times

# Read both datasets
diameters, times = read_times('time-p3m.txt')

# Plotting
plt.plot(diameters, times, marker='o', color='red')

plt.xlabel(r'$r_e/H$')
plt.ylabel('Execution Time (ms)')
#plt.yscale('log')
plt.title('Execution Time vs Cutoff Radius')
plt.grid(True)
#plt.legend()
plt.tight_layout()
plt.show()

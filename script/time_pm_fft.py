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
grid_sizes, times_fftw = read_times('time_pm_fftw.txt')
_, times_kiss = read_times('time_pm_kiss.txt')

ratios = [kiss / fftw for fftw, kiss in zip(times_fftw, times_kiss)]
average_ratio = sum(ratios) / len(ratios)
print(f"Average ratio (KISS / FFTW): {average_ratio:.3f}")

# Plotting
plt.plot(grid_sizes, times_fftw, marker='o', color='red')
plt.plot(grid_sizes, times_kiss, marker='s', color='blue')

plt.xlabel('Grid Size')
plt.ylabel('FFT Time (ms)')
#plt.yscale('log')
plt.title('FFT Time vs Grid Size')
plt.grid(True)
#plt.legend()
plt.tight_layout()
plt.show()

import matplotlib.pyplot as plt
import numpy as np

# File paths
cpu_file = "C:/Projects/ParticleSimulation/pm-timing-cpu.txt"
gpu_file = "C:/Projects/ParticleSimulation/pm-timing-gpu.txt"

def read_timing_data(filepath):
    N_vals = []
    times = []
    with open(filepath, 'r') as file:
        next(file)  # Skip header
        for line in file:
            try:
                n_str, time_str = line.strip().split(',')
                N_vals.append(int(n_str.strip()))
                times.append(float(time_str.strip()))
            except ValueError:
                continue
    return N_vals, times

# Read both files
N_cpu, time_cpu = read_timing_data(cpu_file)
N_gpu, time_gpu = read_timing_data(gpu_file)

# Plot
plt.plot(N_cpu, np.array(time_cpu) / 1000, marker='o', linestyle='-', color='red', label='CPU')
plt.plot(N_gpu, np.array(time_gpu) / 1000, marker='s', linestyle='--', color='green', label='GPU')
plt.title("CPU vs. GPU PM execution time")
plt.xlabel("$N$ (number of particles)")
plt.ylabel("Execution Time (s)")
#plt.yscale("log")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

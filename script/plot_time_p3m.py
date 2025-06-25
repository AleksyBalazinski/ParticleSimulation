import matplotlib.pyplot as plt

def read_data(filename):
    N_values = []
    times = []
    with open(filename, 'r') as file:
        for line in file:
            if line.strip():
                n, t = map(float, line.strip().split())
                N_values.append(n)
                times.append(t)
    return N_values, times

def plot_single(filename, title, color='purple'):
    N, times = read_data(filename)
    N_scaled = [n / 1e4 for n in N]
    times_avg_ms = [t / 100 * 1000 for t in times]  # avg per iteration, convert to ms

    plt.figure()
    plt.plot(N_scaled, times_avg_ms, marker='o', linestyle='-', color=color)
    plt.title(title)
    plt.xlabel(r"$N$ ($\times 10^4$ particles)")
    plt.ylabel("Average Time per Iteration (ms)")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_pm(cpu_file='pm_time_cpu.txt', gpu_file='pm_time_gpu.txt'):
    N_cpu, time_cpu = read_data(cpu_file)
    N_gpu, time_gpu = read_data(gpu_file)

    N_cpu_scaled = [n / 1e4 for n in N_cpu]
    N_gpu_scaled = [n / 1e4 for n in N_gpu]

    time_cpu_avg_ms = [t / 100 * 1000 for t in time_cpu]
    time_gpu_avg_ms = [t / 100 * 1000 for t in time_gpu]

    plt.figure()
    plt.plot(N_cpu_scaled, time_cpu_avg_ms, marker='o', linestyle='-', label='CPU', color='blue')
    plt.plot(N_gpu_scaled, time_gpu_avg_ms, marker='s', linestyle='--', label='GPU', color='green')

    plt.title("PM Simulation Runtime per Iteration")
    plt.xlabel(r"$N$ ($\times 10^4$ particles)")
    plt.ylabel("Average Time per Iteration (ms)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# Example usage:
plot_pm()
plot_single('bh_time.txt', "BH Simulation Runtime per Iteration (CPU)", color='blue')
plot_single('p3m_time.txt', "P3M Simulation Runtime per Iteration (CPU)", color='blue')

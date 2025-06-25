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

import matplotlib.pyplot as plt
import numpy as np

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

def plot_single(filename, title, color='purple', fit_quadratic=False):
    N, times = read_data(filename)
    N_scaled = [n / 1e3 for n in N]
    times_avg_ms = [t / 100 * 1000 for t in times]

    plt.figure()
    plt.plot(N_scaled, times_avg_ms, marker='o', linestyle='-', color=color, label='Measured')

    if fit_quadratic:
        # Fit a quadratic: y = a*x^2 + b*x + c
        coeffs = np.polyfit(N_scaled, times_avg_ms, deg=2)
        poly = np.poly1d(coeffs)
        x_fit = np.linspace(min(N_scaled), max(N_scaled), 200)
        y_fit = poly(x_fit)
        plt.plot(x_fit, y_fit, linestyle='--', color='black', label='Quadratic fit')
        print(f"Quadratic coefficients for {title}: {coeffs}")

    plt.title(title)
    plt.xlabel(r"$N$ ($\times 10^3$ particles)")
    plt.ylabel("Average Time per Iteration (ms)")
    plt.grid(True)
    plt.legend()
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

def plot_bh(standard_file='bh_time_standard.txt', z_order_file='bh_time_z_order.txt'):
    N_cpu, time_cpu = read_data(standard_file)
    N_gpu, time_gpu = read_data(z_order_file)

    N_standard_scaled = [n / 1e4 for n in N_cpu]
    N_z_order_scaled = [n / 1e4 for n in N_gpu]

    time_standard_avg_ms = [t / 100 * 1000 for t in time_cpu]
    time_z_order_avg_ms = [t / 100 * 1000 for t in time_gpu]

    plt.figure()
    plt.plot(N_standard_scaled, time_standard_avg_ms, marker='o', linestyle='-', label='Standard', color='blue')
    plt.plot(N_z_order_scaled, time_z_order_avg_ms, marker='s', linestyle='--', label='Z-ordering', color='green')

    plt.title("PM Simulation Runtime per Iteration")
    plt.xlabel(r"$N$ ($\times 10^4$ particles)")
    plt.ylabel("Average Time per Iteration (ms)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


plot_single('pp_time.txt', "PP Simulation Runtime per Iteration", color='blue', fit_quadratic=True)
import numpy as np
import matplotlib.pyplot as plt

# Constants
DT = 0.04
T = 3
N = 1000
times = np.arange(0, T, DT)

# Physics
def force(theta, g=9.81):
    return -g * np.sin(theta)

def KE(w):
    return 0.5 * w**2

def PE(theta, g=9.81):
    return -g * np.cos(theta)

def area(x, y):
    n = len(x)
    a = 0.0
    for i in range(n):
        j = (i + 1) % n
        a += x[i] * y[j] - y[i] * x[j]
    return abs(a) / 2.0

# Initial square in phase space
x0 = np.concatenate((np.linspace(2, 2, N),
                     np.linspace(2, 1, N),
                     np.linspace(1, 1, N),
                     np.linspace(1, 2, N)))
v0 = np.concatenate((np.linspace(2, 0, N),
                     np.linspace(0, 0, N),
                     np.linspace(0, 2, N),
                     np.linspace(2, 2, N)))

corners_x0 = np.array([2, 2, 1, 1], dtype=float)
corners_v0 = np.array([2, 0, 2, 0], dtype=float)

# Deep copy for both methods
def deep_copy_init():
    return x0.copy(), v0.copy(), corners_x0.copy(), corners_v0.copy()

# Leapfrog
x_lf, v_lf, corners_x_lf, corners_v_lf = deep_copy_init()
areas_lf = []

# half-step kick
v_lf += 0.5 * DT * force(x_lf)
corners_v_lf += 0.5 * DT * force(corners_x_lf)

for t in times:
    x_lf += DT * v_lf
    a = force(x_lf)
    v_lf += DT * a
    v_int_lf = v_lf - 0.5 * DT * a  # interpolated v for area

    corners_x_lf += DT * corners_v_lf
    a_c = force(corners_x_lf)
    corners_v_lf += DT * a_c
    v_int_corners_lf = corners_v_lf - 0.5 * DT * a_c

    areas_lf.append(area(x_lf, v_int_lf))

# Euler
x_eu, v_eu, corners_x_eu, corners_v_eu = deep_copy_init()
areas_eu = []

for t in times:
    areas_eu.append(area(x_eu, v_eu))

    v_new = v_eu + DT * force(x_eu)
    x_eu += DT * v_eu
    v_eu = v_new

    corners_v_new = corners_v_eu + DT * force(corners_x_eu)
    corners_x_eu += DT * corners_v_eu
    corners_v_eu = corners_v_new

# Plotting: Area comparisons
plt.figure(figsize=(8, 5))
plt.plot(times, areas_lf, label="Leapfrog Area", color='blue')
plt.plot(times, areas_eu, label="Euler Area", color='red', linestyle='--')
plt.xlabel("Time")
plt.ylabel("Area in Phase Space")
plt.title("Comparison of Area Conservation: Leapfrog vs Euler")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

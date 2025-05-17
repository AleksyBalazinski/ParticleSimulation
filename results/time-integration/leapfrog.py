import numpy as np
import matplotlib.pyplot as plt

DT = 0.04
T = 3
N = 1000

def force(theta, g=9.81):
    return -g * np.sin(theta)

def KE(w):
    return 0.5 * w**2

def PE(theta, g=9.81):
    return -g * np.cos(theta)

def area(x, y):
    n = len(x)
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += x[i] * y[j]
        area -= y[i] * x[j]
    return abs(area) / 2.0

x0 = np.concatenate((
    np.linspace(2, 2, num=N),
    np.linspace(2, 1, num=N),
    np.linspace(1, 1, num=N),
    np.linspace(1, 2, num=N)))
v0 = np.concatenate((
    np.linspace(2, 0, num=N),
    np.linspace(0, 0, num=N),
    np.linspace(0, 2, num=N),
    np.linspace(2, 2, num=N))) 
x = x0.copy()
v = v0.copy()

corners_x0 = np.array([2, 2, 1, 1])
corners_v0 = np.array([2, 0, 2, 0])
corners_x = corners_x0.copy()
corners_v = corners_v0.copy()

times = []
ke_list = []
pe_list = []
total_energy_list = []
all_positions = [[] for _ in range(len(corners_x))]
all_velocities = [[] for _ in range(len(corners_v))]
areas = []

for i in range(len(corners_x)):
    all_positions[i].append(corners_x[i])
    all_velocities[i].append(corners_v[i])

v = v + 0.5 * DT * force(x) # v(1/2)
corners_v = corners_v + 0.5 * DT * force(corners_x)

for t in np.arange(0, T, DT):
    # leapfrog update
    x = x + DT * v # x(n+1)
    a = force(x)
    v = v + DT * a # v(n+3/2)
    v_int = v - 0.5 * DT * a

    corners_x = corners_x + DT * corners_v
    corners_a = force(corners_x)
    corners_v = corners_v + DT * corners_a
    corners_v_int = corners_v - 0.5 * DT * corners_a

    for i in range(len(corners_x)):
        all_positions[i].append(corners_x[i])
        all_velocities[i].append(corners_v_int[i])

    area_val = area(x, v_int)
    areas.append(area_val)

    ke = KE(corners_v_int[2])
    pe = PE(corners_x[2])
    total_energy = ke + pe

    times.append(t)
    ke_list.append(ke)
    pe_list.append(pe)
    total_energy_list.append(total_energy)


area_max = np.max(areas)
area_min = np.min(areas)
print(f"(max area) - (min area) = {area_max - area_min}")

plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)

plt.plot(times, ke_list, label='Kinetic Energy', linestyle='--')
plt.plot(times, pe_list, label='Potential Energy', linestyle=':')
plt.plot(times, total_energy_list, label='Total Energy', linewidth=2)
plt.xlabel('Time')
plt.ylabel('Energy')
plt.title('Energy Components vs Time')
plt.grid(True)
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(x0, v0)
plt.plot(x, v_int, color='k')

colors = ['c', 'g', 'r', 'm']
for i in range(len(corners_x)):
    plt.plot(all_positions[i], all_velocities[i], color=colors[i], alpha=0.3)
plt.xlabel('Position')
plt.ylabel('Momentum')
plt.title('Phase Space Trajectories')
plt.grid(True)
plt.tight_layout()

plt.figure(figsize=(6, 4))
plt.plot(times, areas, label='Area of Parallelogram', color='purple')
plt.xlabel('Time')
plt.ylabel('Area')
plt.title('Area of Parallelogram in Phase Space vs Time')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.ylim(0.95 * area_min, 1.05 * area_max)

plt.figure(figsize=(8, 4))
for i in range(len(corners_x0)):
    plt.plot(times, all_positions[i][:-1], label=f'Corner {i+1}', color=colors[i])
plt.xlabel('Time')
plt.ylabel('Position (x)')
plt.title('Corner Point Positions vs Time')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
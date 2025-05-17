import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

DT = 0.02
T = 1.0
N = 1000

def force(theta, g=9.81):
    return -g * np.sin(theta)

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

x_history = []
v_history = []

corners_x0 = np.array([2, 2, 1, 1])
corners_v0 = np.array([2, 0, 2, 0])
corners_x = corners_x0.copy()
corners_v = corners_v0.copy()

corner_trajectories_x = [[] for _ in range(4)]
corner_trajectories_v = [[] for _ in range(4)]

for t in np.arange(0, T, DT):
    x_history.append(x.copy())
    v_history.append(v.copy())

    for i in range(4):
        corner_trajectories_x[i].append(corners_x[i])
        corner_trajectories_v[i].append(corners_v[i])

    v_half = v + 0.5 * DT * force(x)
    x = x + DT * v_half
    v = v_half + 0.5 * DT * force(x)

    v_half_corners = corners_v + 0.5 * DT * force(corners_x)
    corners_x = corners_x + DT * v_half_corners
    corners_v = v_half_corners + 0.5 * DT * force(corners_x)

x_history.append(x.copy())
v_history.append(v.copy())
for i in range(4):
    corner_trajectories_x[i].append(corners_x[i])
    corner_trajectories_v[i].append(corners_v[i])

fig, ax = plt.subplots(figsize=(6, 6))
ax.set_xlim(-2.5, 2.5)
ax.set_ylim(-6, 2.5)
ax.set_title("Evolving Quadrilateral with Corner Trajectories")
ax.set_xlabel("x (Position)")
ax.set_ylabel("v (Velocity)")
ax.grid(True)

colors = ['r', 'g', 'm', 'orange']
for i in range(4):
    ax.plot(corner_trajectories_x[i], corner_trajectories_v[i], color=colors[i], lw=1.5, alpha=0.5)

line, = ax.plot([], [], '-', lw=2, color='blue')

def init():
    line.set_data([], [])
    return line,

def update(frame):
    x_vals = x_history[frame]
    v_vals = v_history[frame]

    x_closed = np.append(x_vals, x_vals[0])
    v_closed = np.append(v_vals, v_vals[0])

    line.set_data(x_closed, v_closed)
    return line,

ani = animation.FuncAnimation(
    fig, update, frames=len(x_history),
    init_func=init, blit=True, interval=150, repeat=False
)

plt.tight_layout()
plt.show()

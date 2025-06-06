import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os
import sys
from load_data import *

# Store all trajectories
trajectories = []

frames = load_data(sys.argv[1])
num_particles = len(frames[0])  # Assuming each frame has the same number of particles
trajectories = [np.zeros((len(frames), 2)) for _ in range(num_particles)]

for t, frame in enumerate(frames):
    for i, pos in enumerate(frame):
        trajectories[i][t] = pos[:2]

print(f"x: ({pos_min[0]}, {pos_max[0]}), y: ({pos_min[1]}, {pos_max[1]}), z: ({pos_min[2]}, {pos_max[2]})")
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)

ax.set_xlim(0, 60)
ax.set_ylim(0, 60)
ax.set_xlabel('x (kpc)')
ax.set_ylabel('y (kpc)')
ax.grid(True)

time_step = 1  # 1 Myr per update

def update(frame_num):
    ax.clear()
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 60)
    ax.set_xlabel('x (kpc)')
    ax.set_ylabel('y (kpc)')
    ax.grid(True)
    ax.text(45, 55, f'Time: {frame_num * time_step} Myr', fontsize=10, bbox=dict(facecolor='white', alpha=0.5))
    
    for i in range(num_particles):
        ax.plot(trajectories[i][:frame_num+1, 0], trajectories[i][:frame_num+1, 1], color='blue', alpha=1, linewidth=0.5)
    
    scat = ax.scatter(frames[frame_num][:, 0], frames[frame_num][:, 1], color='red', alpha=1, s=20)
    return scat,

ani = animation.FuncAnimation(fig, update, frames=len(frames), interval=100)
os.makedirs('./animations', exist_ok=True)

def indicate_progress(i, n):
    sys.stdout.write(f'\rSaving frame {i + 1}/{n}')
    sys.stdout.flush()

ani.save('./animations/particles_animation.mp4', writer='ffmpeg', fps=10, progress_callback=lambda i, n: indicate_progress(i, n))

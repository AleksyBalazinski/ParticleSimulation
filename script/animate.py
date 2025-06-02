import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import sys
from load_data import *

frames = load_data(sys.argv[1])
print(f"x: ({pos_min[0]}, {pos_max[0]}), y: ({pos_min[1]}, {pos_max[1]}), z: ({pos_min[2]}, {pos_max[2]})")
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.set_xlim(0, 60)
ax.set_ylim(0, 60)
ax.set_zlim(0, 60)

scat = ax.scatter([], [], [], s=1)

def update(frame):
    ax.clear()
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 60)
    ax.set_zlim(0, 60)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    scat = ax.scatter(frame[:, 0], frame[:, 1], frame[:, 2], color='blue', s=1)
    return scat,

ani = animation.FuncAnimation(fig, update, frames=frames, interval=100)
os.makedirs('./animations', exist_ok=True)

def inicate_progress(i, n):
    sys.stdout.write(f'\rSaving frame {i + 1}/{n}')
    sys.stdout.flush()

ani.save('./animations/particles_animation.mp4', writer='ffmpeg', fps=10, progress_callback=lambda i, n: inicate_progress(i, n))  # Save as MP4
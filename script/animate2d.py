import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import sys
from load_data import *

frames = load_data(sys.argv[1])
print(f"x: ({pos_min[0]}, {pos_max[0]}), y: ({pos_min[1]}, {pos_max[1]}), z: ({pos_min[2]}, {pos_max[2]})")
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)

ax.set_xlim(0, 60)
ax.set_ylim(0, 60)

time_step = 1  # 1 Myr per update
time_elapsed = 0

def update(frame_num):
    global time_elapsed
    ax.clear()
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 60)
    ax.set_xlabel('x (kpc)')
    ax.set_ylabel('y (kpc)')
    ax.text(45, 55, f'Time: {frame_num * time_step} Myr', fontsize=10, bbox=dict(facecolor='white', alpha=0.4))
    scat = ax.scatter(frames[frame_num][:, 0], frames[frame_num][:, 1], color='blue', alpha=0.4, s=0.1)
    return scat,

ani = animation.FuncAnimation(fig, update, frames=len(frames), interval=100)
os.makedirs('./animations', exist_ok=True)

def indicate_progress(i, n):
    sys.stdout.write(f'\rSaving frame {i + 1}/{n}')
    sys.stdout.flush()

ani.save('./animations/particles_animation.mp4', writer='ffmpeg', fps=10, progress_callback=lambda i, n: indicate_progress(i, n))
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os
import sys

def load_data(filename):
    frames = []
    with open(filename, 'r') as f:
        block = []
        for line in f:
            if line.strip():
                block.append(list(map(float, line.split())))
            elif block:  # Blank line after a block
                frames.append(np.array(block))
                block = []
        if block:
            frames.append(np.array(block))
    return frames

frames = load_data(sys.argv[1])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.set_xlim([-110, 110])
ax.set_ylim([-110, 110])
ax.set_zlim([-110, 110])

scat = ax.scatter([], [], [], s=20)

def update(frame):
    ax.clear()
    ax.set_xlim([-110, 110])
    ax.set_ylim([-110, 110])
    ax.set_zlim([-110, 110])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    scat = ax.scatter(frame[:, 0], frame[:, 1], frame[:, 2], color='blue', s=20)
    return scat,

ani = animation.FuncAnimation(fig, update, frames=frames, interval=100)
os.makedirs('./animations', exist_ok=True)

def inicate_progress(i, n):
    sys.stdout.write(f'\rSaving frame {i + 1}/{n}')
    sys.stdout.flush()

ani.save('./animations/particles_animation.mp4', writer='ffmpeg', fps=30, progress_callback=lambda i, n: inicate_progress(i, n))  # Save as MP4
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os
import sys
import math

pos_max = [-math.inf] * 3
pos_min = [math.inf] * 3

def update_limits(position):
    for i in range(3):
        if position[i] > pos_max[i]:
            pos_max[i] = position[i]
        if position[i] < pos_min[i]:
            pos_min[i] = position[i]

def load_data(filename):
    frames = []
    with open(filename, 'r') as f:
        block = []
        for line in f:
            if line.strip():
                position = list(map(float, line.split()))
                block.append(position)
                update_limits(position)
            elif block:  # Blank line after a block
                frames.append(np.array(block))
                block = []
        if block:
            frames.append(np.array(block))
    return frames

frames = load_data(sys.argv[1])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.set_xlim([pos_min[0], pos_max[0]])
ax.set_ylim([pos_min[1], pos_max[1]])
ax.set_zlim([pos_min[2], pos_max[2]])

scat = ax.scatter([], [], [], s=20)

def update(frame):
    ax.clear()
    ax.set_xlim([pos_min[0], pos_max[0]])
    ax.set_ylim([pos_min[1], pos_max[1]])
    ax.set_zlim([pos_min[2], pos_max[2]])
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
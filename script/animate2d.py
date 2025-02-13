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
    print(len(frames))
    return frames

frames = load_data(sys.argv[1])
print(f"x: ({pos_min[0]}, {pos_max[0]}), y: ({pos_min[1]}, {pos_max[1]}), z: ({pos_min[2]}, {pos_max[2]})")
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)

ax.set_xlim(0, 60)
ax.set_ylim(0, 60)

scat = ax.scatter([], [], alpha=0.5, s=0.1)

def update(frame):
    ax.clear()
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 60)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    scat = ax.scatter(frame[:, 0], frame[:, 1], color='blue', alpha=0.5, s=0.1)
    return scat,

ani = animation.FuncAnimation(fig, update, frames=frames, interval=100)
os.makedirs('./animations', exist_ok=True)

def inicate_progress(i, n):
    sys.stdout.write(f'\rSaving frame {i + 1}/{n}')
    sys.stdout.flush()

ani.save('./animations/particles_animation.mp4', writer='ffmpeg', fps=10, progress_callback=lambda i, n: inicate_progress(i, n))  # Save as MP4
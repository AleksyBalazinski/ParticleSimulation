import matplotlib.pyplot as plt
import numpy as np
import math
import sys

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
    print(f"Total frames loaded: {len(frames)}")
    return frames

frames = load_data(sys.argv[1])

last_frame = frames[-1]
x = last_frame[:, 0]
y = last_frame[:, 1]
z = last_frame[:, 2]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(0, 60)
ax.set_ylim(0, 60)
# ax.set_zlim([pos_min[2], pos_max[2]])
ax.set_zlim(0, 60)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

ax.scatter(x, y, z, s=0.1)

plt.show()

import matplotlib.pyplot as plt

plt.figure(figsize=(8, 8))
plt.scatter(x, y, c='blue', marker='o', alpha=0.5, s=0.1)
plt.xlim(0, 60)
plt.ylim(0, 60)  
plt.grid(True) 
plt.show()

plt.scatter(x, z, c='blue', marker='o', alpha=0.5, s=0.1)
plt.xlim(0, 60)
plt.ylim(0, 60)  
plt.grid(True) 
plt.show()

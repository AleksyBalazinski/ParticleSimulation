import matplotlib.pyplot as plt
import numpy as np
import math
import sys

# Initialize limits
pos_max = [-math.inf] * 3
pos_min = [math.inf] * 3

# Function to update position limits
def update_limits(position):
    for i in range(3):
        if position[i] > pos_max[i]:
            pos_max[i] = position[i]
        if position[i] < pos_min[i]:
            pos_min[i] = position[i]

# Function to load data
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

# Load the data from the file specified as a command-line argument
frames = load_data(sys.argv[1])

# Extract the last frame
last_frame = frames[-1]

# Create the plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Set axis limits based on the data
ax.set_xlim([pos_min[0], pos_max[0]])
ax.set_ylim([pos_min[1], pos_max[1]])
ax.set_zlim([pos_min[2], pos_max[2]])

# Set axis labels
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Plot the last frame
ax.scatter(last_frame[:, 0], last_frame[:, 1], last_frame[:, 2], color='blue', s=20)

# Display the plot
plt.show()

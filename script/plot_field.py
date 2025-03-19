import matplotlib.pyplot as plt
import numpy as np
import math
import sys

pos_max = [-math.inf] * 3
pos_min = [math.inf] * 3
m = 1
G = 4.5e-3

def update_limits(position):
    for i in range(3):
        if position[i] > pos_max[i]:
            pos_max[i] = position[i]
        if position[i] < pos_min[i]:
            pos_min[i] = position[i]

def s1_smoothed_g(rs, a):
    gs = []
    for r in rs:
        if r <= a:
            gs.append(-G * m / a**2 * (8 * r / a - 9 * r**2 / a**2 + 2 * r**4 / a**4) )
        else:
            gs.append(-G * m / r**2)

    return gs

def s2_smoothed_g(rs, a):
    gs = []
    for r in rs:
        u = 2 * r / a
        if u <= 1:
            gs.append(-G * m / (35 * a**2) * (224 * u - 224 * u**3 + 70 * u**4 + 48 * u**5 - 21 * u**6))
        elif u <= 2:
            gs.append(-G * m / (35 * a**2) * (12 / u**2 - 224 + 896 * u - 840 * u**2 + 224 * u**3 + 70 * u**4 - 48 * u**5 + 7 * u**6))
        else:
            gs.append(-G * m / r**2)
    
    return gs

def load_data(filename, up_limits = False):
    frames = []
    with open(filename, 'r') as f:
        block = []
        for line in f:
            if line.strip():
                position = list(map(float, line.split()))
                block.append(position)
                if up_limits:
                    update_limits(position)
            elif block:  # Blank line after a block
                frames.append(np.array(block))
                block = []
        if block:
            frames.append(np.array(block))
    print(f"Total frames loaded: {len(frames)}")
    return frames

frames = load_data(sys.argv[1])
field = load_data(sys.argv[2])

last_frame = frames[-1]
x = last_frame[:, 0]
print(min(x))
print(max(x))

last_field = field[-1]
g_x = last_field[:, 0]
g_y = last_field[:, 1]
g_z = last_field[:, 2]

r = (x - 30)
expected_g_x = -4.5e-3 * m/r**2

a = 7.5
smoothed_g_x = s1_smoothed_g(r, a)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(min(x), max(x))
ax.set_ylim(-0.01, 0)


ax.set_xlabel('x (kpc)')
ax.set_ylabel('field value')

ax.plot(x, g_x, label='PM g')
ax.plot(x[5:], expected_g_x[5:], label='expected g')
ax.plot(x, smoothed_g_x, label='smoothed g', linestyle="--")
ax.legend()

plt.show()

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import colormaps
import os
import sys
from load_data import *

if len(sys.argv) != 2:
    print('Missing argument: # of objects')
    sys.exit(1)

option = sys.argv[1]
try:
    number_of_objects = int(option)
    if number_of_objects <= 0:
        raise ValueError()
except ValueError:
    print('Argument must be a positive integer')
    sys.exit(1)

frames = load_data("output/positions.dat")

print(f"x: ({pos_min[0]}, {pos_max[0]}), y: ({pos_min[1]}, {pos_max[1]}), z: ({pos_min[2]}, {pos_max[2]})")

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)

ax.set_xlim(0, number_of_objects * 60)
ax.set_ylim(0, 60)

time_step = 1  # 1 Myr per update

color_cycle = colormaps.get_cmap('tab10').resampled(number_of_objects)

def update(frame_num):
    ax.clear()
    ax.set_xlim(0, number_of_objects * 60)
    ax.set_ylim(0, 60)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('x (kpc)')
    ax.set_ylabel('y (kpc)')
    ax.text(45, 55, f'Time: {frame_num * time_step} Myr', fontsize=10, bbox=dict(facecolor='white', alpha=0.4))

    particles = frames[frame_num]
    total_particles = len(particles)
    chunk_size = total_particles // number_of_objects

    scatters = []
    for i in range(number_of_objects):
        start = i * chunk_size
        end = (i + 1) * chunk_size if i < number_of_objects - 1 else total_particles
        scatters.append(ax.scatter(particles[start:end, 0], particles[start:end, 1], 
                   color=color_cycle(i), alpha=0.2, s=0.4))

    return scatters

ani = animation.FuncAnimation(fig, update, frames=len(frames), interval=100)

os.makedirs('./animations', exist_ok=True)

def indicate_progress(i, n):
    sys.stdout.write(f'\rSaving frame {i + 1}/{n}')
    sys.stdout.flush()

ani.save('./animations/particles_animation.mp4', writer='ffmpeg', fps=10,
         progress_callback=lambda i, n: indicate_progress(i, n))

print("\nAnimation saved at ./animations/particles_animation.mp4")
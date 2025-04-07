import matplotlib.pyplot as plt
import sys
from load_data import *

frames = load_data(sys.argv[1])

last_frame = frames[-1]
x = last_frame[:, 0]
y = last_frame[:, 1]
z = last_frame[:, 2]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(0, 60)
ax.set_ylim(0, 60)
ax.set_zlim(0, 30)

ax.set_xlabel('x (kpc)')
ax.set_ylabel('y (kpc)')
ax.set_zlabel('z (kpc)')

ax.scatter(x, y, z, s=0.1)
plt.show()

plt.figure(figsize=(8, 8))
plt.scatter(x, y, c='blue', marker='o', alpha=0.5, s=0.1)
plt.xlim(0, 60)
plt.ylim(0, 60)
plt.xlabel('x (kpc)')
plt.ylabel('y (kpc)')
plt.grid(True) 
plt.show()

plt.scatter(x, z, c='blue', marker='o', alpha=0.5, s=0.1)
plt.xlim(0, 60)
plt.ylim(0, 60)
plt.xlabel('x (kpc)')
plt.ylabel('z (kpc)')
plt.grid(True) 
plt.show()

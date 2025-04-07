import matplotlib.pyplot as plt
import numpy as np
import sys
from load_data import *

m = 1
G = 4.5e-3

def s1_smoothed_g(rs, a):
    gs = []
    for r in rs:
        if r <= a:
            gs.append(G * m / a**2 * (8 * r / a - 9 * r**2 / a**2 + 2 * r**4 / a**4) )
        else:
            gs.append(G * m / r**2)

    return gs

def s2_smoothed_g(rs, a):
    gs = []
    for r in rs:
        u = 2 * r / a
        if u <= 1:
            gs.append(G * m / (35 * a**2) * (224 * u - 224 * u**3 + 70 * u**4 + 48 * u**5 - 21 * u**6))
        elif u <= 2:
            gs.append(G * m / (35 * a**2) * (12 / u**2 - 224 + 896 * u - 840 * u**2 + 224 * u**3 + 70 * u**4 - 48 * u**5 + 7 * u**6))
        else:
            gs.append(G * m / r**2)
    
    return gs

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
expected_g_x = G* m/r**2
H = 1.875
a = 2.5 * H
smoothed_g_x = s1_smoothed_g(r, a)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(0, 10)
ax.set_ylim(0, 0.15)


ax.set_xlabel('$r/H$')
ax.set_ylabel('$F / (Gm_1m_2)$')

ax.plot(np.array(r)/H, -np.array(g_x)/G, label='PM', linestyle="-", color="blue", linewidth=2)
ax.plot(np.array(r[5:]) / H, np.array(expected_g_x[5:])/G, label='Inverse square', linestyle=":", color="red", linewidth=2)
ax.plot(np.array(r) / H, np.array(smoothed_g_x)/G, label='Reference', linestyle="--", color="green", linewidth=2)
ax.legend()

plt.show()

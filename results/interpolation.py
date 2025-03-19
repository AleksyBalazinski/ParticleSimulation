import matplotlib.pyplot as plt
import numpy as np

x_squared_range = np.linspace(0, np.pi ** 2, 1001)
delta_squared = x_squared_range[1] - x_squared_range[0]
sin_x_table = np.sin(np.sqrt(x_squared_range)) # precomputed 

xs = np.linspace(0, np.pi, endpoint=False, num=100)
xs_squared = np.square(xs) # hypothetical input
sin_interpolated = []
for x_squared in xs_squared:
    ksi = x_squared / delta_squared
    t = int(ksi)
    sin_x = sin_x_table[t] + (ksi - t) * (sin_x_table[t+1] - sin_x_table[t])
    sin_interpolated.append(sin_x)


plt.scatter(xs, sin_interpolated, s=10)
plt.plot(xs, np.sin(xs), linestyle="--")
plt.show()
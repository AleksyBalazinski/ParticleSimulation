import matplotlib.pyplot as plt
import numpy as np

def sample_from_disk(rd, n):
    phi = np.random.random(n) * 2 * np.pi
    cdf = np.random.random(n)
    r = rd * np.sqrt(cdf)

    x = r * np.cos(phi)
    y = r * np.sin(phi)

    return x, y

def get_count_in_box(x, y, xc, yc, dx, dy):
    cnt = 0
    for p in zip(x, y):
        if p[0] >= xc - dx/2 and p[0] <= xc + dx/2 \
            and p[1] >= yc - dy/2 and p[1] <= yc + dy/2:
            cnt += 1

    return cnt

def get_counts_along_x(x, y, rd, n):
    counts = np.zeros(n)
    xcs = np.zeros(n)
    dx = rd / n
    for i in range(n):
        xc = i * dx + dx / 2
        xcs[i] = xc
        counts[i] = get_count_in_box(x, y, xc, 0, .1, .1)

    return xcs, counts

x, y = sample_from_disk(1, int(1e6))

plt.figure(figsize=(8, 8))
plt.scatter(x, y, c='blue', marker='o', alpha=0.7, s=1)

# Labels and title
plt.xlabel('X Axis')
plt.ylabel('Y Axis')
plt.title('2D Scatter Plot')

# Show plot
plt.grid(True)  # Adds grid for better visibility
plt.show()

centers, counts = get_counts_along_x(x, y, 1, 50)
plt.plot(centers, counts)
plt.show()
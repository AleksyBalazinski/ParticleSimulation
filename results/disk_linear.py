import matplotlib.pyplot as plt
import numpy as np

def my_func(r, rd, f):
    return 2 * r**3 - 3 * rd * r**2 + f * rd**3

def my_func_der(r, rd):
    return 6 * r**2 - 6 * rd * r

def newton(f, df, x0, max_iter, eps):
    i = 0
    err = np.inf
    x = x0
    while err > eps and i < max_iter:
        i += 1
        x_prev = x
        x = x - f(x) / df(x)
        err = np.abs(x - x_prev)
    
    return x

def sample_from_linear(rd, n):
    nums = np.zeros(n)
    for i in range(n):
        cdf = np.random.random()
        f = lambda r: my_func(r, rd, cdf)
        df = lambda r: my_func_der(r, rd)
        nums[i] = newton(f, df, rd / 2, 100, 0.01)

    return nums

def sample_from_disk(rd, n):
    phi = np.random.random(n) * 2 * np.pi
    r = sample_from_linear(rd, n)

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


x, y = sample_from_disk(1, int(1e5))

plt.figure(figsize=(8, 8))
plt.scatter(x, y, c='blue', marker='o', alpha=0.7, s=.1)

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

import matplotlib.pyplot as plt
import numpy as np

def myFunc(r, rD, f):
    return 3 * r**4 - 4 * rD * r**3 + f * rD**4 

def myFuncDer(r, rD, f):
    return 12 * r**3 - 12 * rD * r**2

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
        f = np.random.random()
        func = lambda r: myFunc(r, rd, f)
        funcDer = lambda r: myFuncDer(r, rd, f)
        nums[i] = newton(func, funcDer, rd / 2, 100, 0.01)
        
    return nums

def sample_from_spherical(rd, n):
    phi = np.random.random(n) * 2 * np.pi
    costheta = np.random.random(n) * 2 - 1
    theta = np.arccos(costheta)
    r = sample_from_linear(rd, n)

    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    return x, y, z

def get_count_in_box(x, y, z, xc, yc, zc, dx, dy, dz):
    cnt = 0
    for p in zip(x, y, z):
        if p[0] >= xc - dx/2 and p[0] <= xc + dx/2 \
            and p[1] >= yc - dy/2 and p[1] <= yc + dy/2 \
            and p[2] >= zc - dz/2 and p[2] <= zc + dz/2:
            cnt += 1

    return cnt

def get_counts_along_x(x, y, z, rd, n):
    counts = np.zeros(n)
    xcs = np.zeros(n)
    dx = rd / n
    for i in range(n):
        xc = i * dx + dx / 2
        xcs[i] = xc
        counts[i] = get_count_in_box(x, y, z, xc, 0, 0, .1, .1, .1)

    return xcs, counts


x, y, z = sample_from_spherical(2, int(1e4))

plt.figure(figsize=(8, 8))
plt.scatter(x, z, c='blue', marker='o', alpha=0.7, s=1)

plt.xlabel('X Axis')
plt.ylabel('Y Axis')
plt.title('2D Scatter Plot')

plt.grid(True)
plt.show()

centers, counts = get_counts_along_x(x, y, z, 1, 50)
plt.plot(centers, counts)
plt.show()
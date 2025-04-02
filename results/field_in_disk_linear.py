import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import mpmath

def gx(x, u, phi):
    eps = 0.001
    denom = (x**2 + u**2 + 2 * x * u * np.cos(phi) + eps**2)**(3/2)
    return u * (x + u * np.cos(phi)) / denom

def gx_linear(x, u, phi):
    eps = 0.001
    denom = (x**2 + u**2 + 2 * x * u * np.cos(phi) + eps**2)**(3/2)
    return u * (1 - u) * (x + u * np.cos(phi)) / denom

def gx_paper(x):
    n = 4 * x / (1 + x)**2
    k = np.sqrt(n)
    eps = 0
    term1 = -((1 + x**2) / (x * (1 + x) + eps)) * mpmath.ellipk(n)
    term2 = (1 + x) * (2 + x) / (2 * x + eps) * mpmath.ellipe(n)
    term3 = -(1 - x)**2 / (2 * (1 + x)) * mpmath.ellippi(n, n)
    return -2 * (term1 + term2 + term3)

n = 50
x_space = np.linspace(0.01, 1.5, n)
g_space = np.zeros(n)
for i, x in enumerate(x_space):
    f = lambda u, phi: gx(x, u, phi)
    g_space[i] =  integrate.dblquad(f, 0, 2 * np.pi, 0, 1)[0]

g_space_linear = np.zeros(n)
for i, x in enumerate(x_space):
    f = lambda u, phi: gx_linear(x, u, phi)
    g_space_linear[i] =  integrate.dblquad(f, 0, 2 * np.pi, 0, 1)[0]

g_space_paper = np.zeros(n)
for i, x in enumerate(x_space):
    g_space_paper[i] =  gx_paper(x)

k = 2.5
h = .66
a = -k / h**2
g_space_approx = a * (x_space - h)**2 + k  

plt.plot(x_space, g_space_paper, label="Constant Density", linestyle="--", color="blue", linewidth=2)
plt.plot(x_space, g_space_linear, label="Radially Decreasing Density", linestyle="-", color="red", linewidth=2)
plt.plot(x_space, g_space_approx, label="Radially Decreasing Density (Approx.)", linestyle=":", color="green", linewidth=2)

plt.grid(True, linestyle="--", alpha=0.5)

plt.xlabel(r"$r/R$")
plt.ylabel(r"$g_r/(G\sigma_0)$")

plt.legend(loc="best", frameon=True)

plt.show()



import numpy as np
import matplotlib.pyplot as plt

Norm = 512 / (7 * np.pi)

def F(x): 
    sqrt_term = np.sqrt(1 - x**2)
    arctan_term = np.atan(x / sqrt_term)
    numerator = (
        x * sqrt_term * (-105 + 1210 * x**2 - 2104 * x**4 + 1488 * x**6 - 384 * x**8)
        + 105 * arctan_term
    )
    return numerator / 3840 * Norm

def p(x): 
    return x**2 * (1 - x**2)**(7/2) * Norm

def inverse_transform_newton(r, tol=1e-8, max_iter=100):
    x = 0.5 
    for _ in range(max_iter):
        fx = F(x) - r
        dfx = p(x)
        if abs(fx) < tol:
            return x
        if dfx == 0:
            break
        x = x - fx / dfx
    return x

np.random.seed(42)
N = 10000
uniform_samples = np.random.uniform(0, 1, N)
x_samples = np.array([inverse_transform_newton(r) for r in uniform_samples])

plt.hist(x_samples, bins=25, density=True, alpha=0.6, label='Sampled')

x = np.linspace(0, 1, 500)
pdf_vals = p(x)
plt.plot(x, pdf_vals, 'r-', label='True PDF')

plt.title('Inverse Transform Sampling')
plt.xlabel('x')
plt.ylabel('Density')
plt.legend()
plt.grid(True)
plt.show()

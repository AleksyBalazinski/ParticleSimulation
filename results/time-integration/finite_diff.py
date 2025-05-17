import numpy as np
import matplotlib.pyplot as plt

H_values = np.logspace(-6, 0, 50)
errors_basic = []
errors_improved = []

a = 4 / 3

for H in H_values:
    x = np.arange(-2 * np.pi, 2 * np.pi, H)
    
    basic_derivative = (np.sin(x + H) - np.sin(x - H)) / (2 * H)
    
    improved_derivative = (
        a * (np.sin(x + H) - np.sin(x - H)) / (2 * H) +
        (1 - a) * (np.sin(x + 2*H) - np.sin(x - 2*H)) / (4 * H)
    )
    
    exact_derivative = np.cos(x)
    
    errors_basic.append(np.max(np.abs(basic_derivative - exact_derivative)))
    errors_improved.append(np.max(np.abs(improved_derivative - exact_derivative)))

H = 0.1
a_values = np.linspace(-1, 2, 100)
errors_vs_a = []

x = np.arange(-2 * np.pi, 2 * np.pi, H)
exact_derivative = np.cos(x)

for a in a_values:
    improved_derivative = (
        a * (np.sin(x + H) - np.sin(x - H)) / (2 * H) +
        (1 - a) * (np.sin(x + 2*H) - np.sin(x - 2*H)) / (4 * H)
    )
    error = np.max(np.abs(improved_derivative - exact_derivative))
    errors_vs_a.append(error)

fig, axs = plt.subplots(1, 2, figsize=(14, 5))

axs[0].loglog(H_values, errors_basic, '-', label='Basic (2nd order)')
axs[0].loglog(H_values, errors_improved, '-', label='Improved (4th order)', color='orange')
axs[0].set_title('Max Absolute Error')
axs[0].set_xlabel('H (step size)')
axs[0].set_ylabel('Error')
axs[0].grid(True, ls='--')
axs[0].legend()

axs[1].plot(a_values, errors_vs_a, '-')
axs[1].set_title(r'Error of Improved Finite Difference vs Parameter $\alpha$')
axs[1].set_xlabel(r'$\alpha$')
axs[1].set_ylabel('Max Absolute Error')
axs[1].set_yscale('log')
axs[1].grid(True, ls='--')

plt.tight_layout()
plt.show()

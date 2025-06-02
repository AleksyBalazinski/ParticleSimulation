import numpy as np
import matplotlib.pyplot as plt

# Load data from file (adjust delimiter if needed: ',' or '\t')
data = np.loadtxt('time-p3m.txt', delimiter=',')
x = data[:, 0]  # Time values
y = data[:, 1]  # Corresponding measurements

# Fit a polynomial of specified degree
degree = 2  # Change as needed
coefficients = np.polyfit(x, y, degree)
polynomial = np.poly1d(coefficients)

# Generate fitted curve
x_fit = np.linspace(min(x), max(x), 500)
y_fit = polynomial(x_fit)

# Plotting
plt.scatter(x, y, color='red', label='Data points')
plt.plot(x_fit, y_fit, label=f'{degree}° Polynomial fit')
plt.xlabel('Time')
plt.ylabel('Value')
plt.title('Polynomial Fit to Data in time-p3m.txt')
plt.legend()
plt.grid(True)
plt.show()

# Print polynomial equation
print("Fitted polynomial coefficients:", coefficients)
print("Polynomial equation:\n", polynomial)

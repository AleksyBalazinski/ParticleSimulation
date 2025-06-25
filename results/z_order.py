import matplotlib.pyplot as plt

def morton_curve_interleaved_yx(n=2):
    """
    Generate Morton curve points with bit interleaving order:
    bit pair = (Y bit, then X bit), meaning X bit is the LSB in each pair.
    """
    size = 2 ** n
    points = []

    for i in range(size * size):
        x = y = 0
        for j in range(n):
            # Now reverse the order of bits in the pair: Y bit first, X bit second
            y |= ((i >> (2 * j + 1)) & 1) << j
            x |= ((i >> (2 * j)) & 1) << j
        points.append((x, y, i))  # include Morton code
    return points

def plot_colored_binary(ax, x, y, code, n=2, dx=0.1):
    """Plot each bit of code at (x,y) with Y bits black and X bits red,
    with bits ordered as Y (even indices), X (odd indices)."""
    bits = format(code, f'0{2*n}b')
    for i, bit in enumerate(bits):
        # Since bit pairs are (Y,X), even bits = Y (black), odd bits = X (red)
        color = 'black' if i % 2 == 0 else 'red'
        ax.text(x + i*dx, y, bit, fontsize=12, color=color, va='bottom', ha='left')

# Parameters
n = 2
size = 2 ** n
points = morton_curve_interleaved_yx(n=n)

fig, ax = plt.subplots(figsize=(6, 6))

# Plot Z-order curve lines
for i in range(len(points) - 1):
    x0, y0, _ = points[i]
    x1, y1, _ = points[i + 1]
    # Reverse y-axis by plotting with flipped y: y -> size - 1 - y
    ax.plot([x0, x1], [size - 1 - y0, size - 1 - y1], color='blue')

# Plot points and colored Morton codes (with reversed y)
for x, y, code in points:
    ax.plot(x, size - 1 - y, 'bo', markersize=4)
    plot_colored_binary(ax, x + 0.05, size - 1 - y + 0.05, code, n=n)

# Label axes with colored binary labels
for i in range(size):
    x_label = format(i, f'0{n}b')
    y_label = format(i, f'0{n}b')
    # X axis labels in red (since now X bits are odd bits, let's keep this consistent)
    ax.text(i, -0.6, x_label, ha='center', fontsize=12, color='red')
    # Y axis labels in black, reversed order on axis
    ax.text(-0.6, size - 1 - i, y_label, va='center', fontsize=12, color='black')

ax.set_xlim(-1, size)
ax.set_ylim(-1, size - 0.5)
ax.set_xticks(range(size))
ax.set_yticks(range(size))

# Reverse the y-axis labels only (not the axis direction)
yticks = range(size)
ax.set_yticks(yticks)
ax.set_yticklabels(list(yticks)[::-1])

ax.grid(True, linestyle='--', linewidth=0.5)
ax.set_aspect('equal', adjustable='box')
plt.tight_layout()
plt.show()


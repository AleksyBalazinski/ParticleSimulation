import matplotlib.pyplot as plt

# Thread counts
thread_counts = [1, 4, 8, 12]
colors = ['blue', 'green', 'orange', 'red']

plt.figure()

# Loop through each file
for threads, color in zip(thread_counts, colors):
    filename = f'p3m_{threads}_thread_time.txt'
    particles = []
    time_ms = []
    
    # Read data
    with open(filename, 'r') as f:
        for line in f:
            if line.strip():  # Skip empty lines
                n, t = line.strip().split()
                particles.append(int(n) / 1e4)
                time_ms.append(float(t) / 25)  # Convert seconds to milliseconds

    # Plot
    plt.plot(particles, time_ms, label=(f'{threads} Thread' + ('s' if threads > 1 else '')), color=color, marker='o')

# Labeling
plt.title('P3M Running Time vs Number of Particles')
plt.xlabel(r'Number of Particles ($\times 10^4$)')
plt.ylabel('Running Time (s)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

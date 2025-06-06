import numpy as np
import matplotlib.pyplot as plt
import sys
from load_data import *

if __name__ == "__main__":
    frames1 = load_data(sys.argv[1])
    frames2 = load_data(sys.argv[2])

    if len(frames1) != len(frames2):
        raise ValueError("The files don't have equal length")
    if len(frames1[0]) != len(frames2[0]):
        raise ValueError("Numbers of particles are not the same in both files")
    
    differences = np.zeros(len(frames1))
    i = 0
    for f1, f2 in zip(frames1, frames2):
        for p1, p2 in zip(f1, f2):
            differences[i] += (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2
        i += 1

    plt.plot(differences)
    plt.title('Trajectory difference')
    plt.show()
    


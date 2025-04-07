import numpy as np
import math
import struct

pos_max = [-math.inf] * 3
pos_min = [math.inf] * 3

def update_limits(position):
    for i in range(3):
        if position[i] > pos_max[i]:
            pos_max[i] = position[i]
        if position[i] < pos_min[i]:
            pos_min[i] = position[i]

def load_data(filename):
    frames = []
    with open(filename, 'rb') as f:
        header = f.read(8)
        particles_cnt, frames_cnt = struct.unpack("<ii", header)
        total_vectors = particles_cnt * frames_cnt
        total_floats = total_vectors * 3
        float_data = struct.unpack(f"<{total_floats}f", f.read(total_floats * 4))
        block = []
        cur_vector = []
        for x in float_data:
            cur_vector.append(x)
            if len(cur_vector) == 3:
                block.append(cur_vector)
                update_limits(cur_vector)
                cur_vector = []
            if len(block) == particles_cnt:
                frames.append(np.array(block))
                block = []
            if len(frames) == frames_cnt:
                break
    print(len(frames))
    return frames
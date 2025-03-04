#pragma once

#include <vector>
#include "vec3.h"

std::string ppMethod(std::vector<Vec3>& state,
                     std::vector<float>& masses,
                     const int simLength,
                     const float stepSize,
                     const float G,
                     const char* positionsPath = "output.txt",
                     const char* energyPath = "energy.txt",
                     const char* momentumPath = "momentum.txt");
#pragma once

#include <vector>
#include "vec3.h"

std::string ppMethod(std::vector<Vec3>& state,
                     std::vector<double>& masses,
                     const int simLength,
                     const double stepSize,
                     const double G,
                     const char* positionsPath = "output.txt",
                     const char* energyPath = "energy.txt",
                     const char* momentumPath = "momentum.txt");
#pragma once

#include <vector>
#include "vec3.h"

std::string ppMethodRK4(std::vector<Vec3>& state,
                        std::vector<float>& masses,
                        const int simLength,
                        const float stepSize,
                        const float G,
                        const char* positionsPath = "output.txt",
                        const char* energyPath = "energy.txt",
                        const char* momentumPath = "momentum.txt");

std::string ppMethodLeapfrog(const std::vector<Vec3>& state,
                             const std::vector<float>& masses,
                             const int simLength,
                             const float stepSize,
                             const float G,
                             const char* positionsPath = "output.txt",
                             const char* energyPath = "energy.txt",
                             const char* momentumPath = "momentum.txt",
                             const char* angularMomentumPath = "angular_momentum.txt");
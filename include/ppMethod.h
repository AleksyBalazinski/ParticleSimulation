#pragma once

#include <string>
#include <vector>
#include "RK4Stepper.h"
#include "stateRecorder.h"
#include "vec3.h"

std::string ppMethod(std::vector<Vec3>& state,
                     std::vector<double>& masses,
                     const double simLengthSeconds = 10.0,
                     const double stepSize = 0.001,
                     const double G = 1,
                     const int frameRate = 30,
                     const char* positionsPath = "output.txt",
                     const char* energyPath = "energy.txt",
                     const char* momentumPath = "momentum.txt");
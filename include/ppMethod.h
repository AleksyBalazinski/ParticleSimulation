#pragma once

#include <vector>
#include "stateRecorder.h"
#include "vec3.h"

std::string ppMethodRK4(std::vector<Vec3>& state,
                        std::vector<float>& masses,
                        const int simLength,
                        const float stepSize,
                        const float G,
                        StateRecorder& stateRecorder);

std::string ppMethodLeapfrog(const std::vector<Vec3>& state,
                             const std::vector<float>& masses,
                             const int simLength,
                             const float stepSize,
                             const float G,
                             const float softeningLength,
                             StateRecorder& stateRecorder,
                             bool recordField = false);
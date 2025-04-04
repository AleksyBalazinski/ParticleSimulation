#pragma once

#include <vector>
#include "triple.cuh"
#include "utils.cuh"
#include "vec3.cuh"

void pmMethod(std::vector<Vec3>& state,
              const std::vector<float>& masses,
              Triple<float> effectiveBoxSize,
              float H,
              float DT,
              SphRadDecrFieldParams params,
              float G,
              int simLength);
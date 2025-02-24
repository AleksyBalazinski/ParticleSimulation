#pragma once

#include <vector>
#include "ext_field_bind.cuh"
#include "vec3.cuh"

void pmMethod(std::vector<Vec3>& state,
              const std::vector<float>& masses,
              float effectiveBoxSize,
              float H,
              float DT,
              float G,
              int simLength);
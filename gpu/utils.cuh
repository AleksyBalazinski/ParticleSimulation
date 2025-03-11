#pragma once

#include <functional>
#include <vector>
#include "vec3.cuh"

#define PI 3.1415926535897932f

__host__ float newton(std::function<float(float)> f,
                      std::function<float(float)> df,
                      float guess,
                      int maxIter,
                      float tolerance);
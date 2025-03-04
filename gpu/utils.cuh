#pragma once

#include <functional>
#include <vector>
#include "vec3.cuh"

#define PI 3.1415926535897932f

__device__ Vec3 externalFieldBulge(Vec3 pos, Vec3 bulge, float rb, float mb, float G);

__host__ Vec3 externalFieldBulgeHost(Vec3 pos, Vec3 bulge, float rb, float mb, float G);

__host__ float newton(std::function<float(float)> f,
                      std::function<float(float)> df,
                      float guess,
                      int maxIter,
                      float tolerance);
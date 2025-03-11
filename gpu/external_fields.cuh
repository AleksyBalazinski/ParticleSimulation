#pragma once
#include <functional>
#include "vec3.cuh"

struct SphRadDecrFieldParams {
  __host__ __device__ SphRadDecrFieldParams(Vec3 center, float r, float m);
  Vec3 center;
  float r;
  float m;
};

__device__ Vec3 sphRadDecrField(Vec3 pos, SphRadDecrFieldParams params, float G);

__host__ Vec3 sphRadDecrFieldHost(Vec3 pos, SphRadDecrFieldParams params, float G);

__host__ float newton(std::function<float(float)> f,
                      std::function<float(float)> df,
                      float guess,
                      int maxIter,
                      float tolerance);
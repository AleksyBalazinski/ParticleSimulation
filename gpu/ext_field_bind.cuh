#pragma once

#include "vec3.cuh"

class ExtFieldBind {
 private:
  Vec3 center;
  float rb;
  float mb;
  float G;
  Vec3 (*field)(Vec3, Vec3, float, float, float);

 public:
  ExtFieldBind(Vec3 center,
               float rb,
               float mb,
               float G,
               Vec3 (*field)(Vec3, Vec3, float, float, float));
  __host__ __device__ Vec3 operator()(Vec3 pos);
};
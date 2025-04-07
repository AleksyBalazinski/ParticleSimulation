#pragma once

#include "vec3.h"

inline __device__ Vec3 makeVec3(float x, float y, float z) {
  Vec3 v;
  v.x = x;
  v.y = y;
  v.z = z;

  return v;
}

inline __device__ Vec3 add(Vec3 a, Vec3 b) {
  return makeVec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline __device__ Vec3& increment(Vec3& target, Vec3 a) {
  target.x += a.x;
  target.y += a.y;
  target.z += a.z;

  return target;
}

template <typename... Ts>
__device__ Vec3 add(Vec3 head, Ts... tail) {
  return add(head, add(tail...));
}

inline __device__ Vec3 sub(Vec3 a, Vec3 b) {
  return makeVec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline __device__ Vec3 mul(float s, Vec3 a) {
  return makeVec3(s * a.x, s * a.y, s * a.z);
}

inline __device__ Vec3 div(Vec3 a, float s) {
  return makeVec3(a.x / s, a.y / s, a.z / s);
}
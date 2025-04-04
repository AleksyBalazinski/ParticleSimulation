#pragma once

template <typename T>
struct Triple {
  __device__ __host__ Triple(T x, T y, T z) : x(x), y(y), z(z) {}
  T x;
  T y;
  T z;
};
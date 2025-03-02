#pragma once

#include <string>

class Vec3 {
 public:
  float x;
  float y;
  float z;

  __host__ __device__ Vec3(float x, float y, float z);

  __host__ __device__ Vec3();

  __host__ __device__ Vec3& operator+=(const Vec3 other);

  __host__ __device__ float getMagnitude() const;

  __host__ __device__ float getMagnitudeSquared() const;

  __host__ char* Vec3::toString(char* singleBuf,
                                std::size_t singleBufSize,
                                char* vecBuf,
                                std::size_t vecBufSize) const;

  __host__ __device__ Vec3(const float3& from);

  __host__ __device__ operator float3();
};

__host__ __device__ Vec3 operator+(const Vec3& a, const Vec3& b);

__host__ __device__ Vec3 operator-(const Vec3& a, const Vec3& b);

__host__ __device__ Vec3 operator*(float s, const Vec3& a);

__host__ __device__ Vec3 operator/(const Vec3& a, float s);
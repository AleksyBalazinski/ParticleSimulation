#include <cmath>
#include <cstdio>
#include <cstring>
#include "vec3.cuh"

__host__ Vec3 operator+(const Vec3& a, const Vec3& b) {
  return Vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

__host__ Vec3 operator-(const Vec3& a, const Vec3& b) {
  return Vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

__host__ Vec3 operator*(float s, const Vec3& a) {
  return Vec3(s * a.x, s * a.y, s * a.z);
}

__host__ Vec3 operator/(const Vec3& a, float s) {
  return Vec3(a.x / s, a.y / s, a.z / s);
}

__host__ __device__ Vec3::Vec3(float x, float y, float z) : x(x), y(y), z(z) {}

__host__ __device__ Vec3::Vec3() : x(0), y(0), z(0) {}

__host__ Vec3& Vec3::operator+=(const Vec3 other) {
  x += other.x;
  y += other.y;
  z += other.z;

  return *this;
}

__host__ float Vec3::getMagnitude() const {
  return std::sqrt(x * x + y * y + z * z);
}

__host__ float Vec3::getMagnitudeSquared() const {
  return x * x + y * y + z * z;
}

__host__ char* Vec3::toString(char* singleBuf,
                              std::size_t singleBufSize,
                              char* vecBuf,
                              std::size_t vecBufSize) const {
  std::memset(vecBuf, 0, vecBufSize);
  int cnt1 = std::snprintf(singleBuf, singleBufSize, "%f", x);
  std::memcpy(vecBuf, singleBuf, cnt1);
  vecBuf[cnt1] = ' ';

  int cnt2 = std::snprintf(singleBuf, singleBufSize, "%f", y);
  std::memcpy(vecBuf + cnt1 + 1, singleBuf, cnt2);
  vecBuf[cnt1 + 1 + cnt2] = ' ';

  int cnt3 = std::snprintf(singleBuf, singleBufSize, "%f", z);
  std::memcpy(vecBuf + cnt1 + cnt2 + 2, singleBuf, cnt3);

  return vecBuf;
}

__host__ __device__ Vec3::Vec3(const float3& from) : x(from.x), y(from.y), z(from.z) {}

__host__ __device__ Vec3::operator float3() {
  return make_float3(x, y, z);
};

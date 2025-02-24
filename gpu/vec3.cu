#include <cmath>
#include "vec3.cuh"

__host__ __device__ Vec3 operator+(const Vec3& a, const Vec3& b) {
  return Vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

__host__ __device__ Vec3 operator-(const Vec3& a, const Vec3& b) {
  return Vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

__host__ __device__ Vec3 operator*(float s, const Vec3& a) {
  return Vec3(s * a.x, s * a.y, s * a.z);
}

__host__ __device__ Vec3 operator/(const Vec3& a, float s) {
  return Vec3(a.x / s, a.y / s, a.z / s);
}

__host__ __device__ Vec3::Vec3(float x, float y, float z) : x(x), y(y), z(z) {}

__host__ __device__ Vec3::Vec3() : x(0), y(0), z(0) {}

__host__ __device__ Vec3& Vec3::operator+=(const Vec3 other) {
  x += other.x;
  y += other.y;
  z += other.z;

  return *this;
}

__host__ __device__ float Vec3::getMagnitude() const {
  return std::sqrt(x * x + y * y + z * z);
}

__host__ __device__ float Vec3::getMagnitudeSquared() const {
  return x * x + y * y + z * z;
}

__host__ std::string Vec3::toString() const {
  return std::to_string(x) + ' ' + std::to_string(y) + ' ' + std::to_string(z);
}

__host__ __device__ Vec3::Vec3(const float3& from) : x(from.x), y(from.y), z(from.z) {}

__host__ __device__ Vec3::operator float3() {
  return make_float3(x, y, z);
};

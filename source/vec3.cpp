#include "vec3.h"
#include <cmath>
#include <cstdio>
#include <cstring>

Vec3 operator+(const Vec3& a, const Vec3& b) {
  return Vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

Vec3 operator-(const Vec3& a, const Vec3& b) {
  return Vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

Vec3 operator*(float s, const Vec3& a) {
  return Vec3(s * a.x, s * a.y, s * a.z);
}

Vec3 operator/(const Vec3& a, float s) {
  return Vec3(a.x / s, a.y / s, a.z / s);
}

Vec3 Vec3::create(float x, float y, float z) {
  Vec3 v;
  v.x = x;
  v.y = y;
  v.z = z;

  return v;
}

Vec3& Vec3::operator+=(const Vec3 other) {
  x += other.x;
  y += other.y;
  z += other.z;

  return *this;
}

Vec3& Vec3::operator/=(float s) {
  x /= s;
  y /= s;
  z /= s;
  return *this;
}

Vec3 Vec3::cross(const Vec3 other) const {
  return Vec3(y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x);
}

float Vec3::getMagnitude() const {
  return std::sqrt(x * x + y * y + z * z);
}

float Vec3::getMagnitudeSquared() const {
  return x * x + y * y + z * z;
}

char* Vec3::toString(char* singleBuf,
                     std::size_t singleBufSize,
                     char* vecBuf,
                     std::size_t vecBufSize) const {
  std::memset(vecBuf, 0, vecBufSize);
  int cnt1 = std::snprintf(singleBuf, singleBufSize, "%.6e", x);
  std::memcpy(vecBuf, singleBuf, cnt1);
  vecBuf[cnt1] = ' ';

  int cnt2 = std::snprintf(singleBuf, singleBufSize, "%.6e", y);
  std::memcpy(vecBuf + cnt1 + 1, singleBuf, cnt2);
  vecBuf[cnt1 + 1 + cnt2] = ' ';

  int cnt3 = std::snprintf(singleBuf, singleBufSize, "%.6e", z);
  std::memcpy(vecBuf + cnt1 + cnt2 + 2, singleBuf, cnt3);

  return vecBuf;
}

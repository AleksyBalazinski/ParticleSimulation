#pragma once

#include <string>

struct Vec3 {
  float x{};
  float y{};
  float z{};

  static Vec3 create(float x, float y, float z);

  Vec3& operator+=(const Vec3 other);

  Vec3& operator/=(float s);

  Vec3 cross(const Vec3 other) const;

  float getMagnitude() const;

  float getMagnitudeSquared() const;

  char* toString(char* singleBuf, size_t singleBufSize, char* vecBuf, size_t vecBufSize) const;

  static Vec3 zero() {
    Vec3 z;
    z.x = 0;
    z.y = 0;
    z.z = 0;

    return z;
  }
};

Vec3 operator+(const Vec3& a, const Vec3& b);

Vec3 operator-(const Vec3& a, const Vec3& b);

Vec3 operator*(float s, const Vec3& a);

Vec3 operator/(const Vec3& a, float s);
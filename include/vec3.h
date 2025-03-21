#pragma once

#include <string>

class Vec3 {
 public:
  float x;
  float y;
  float z;

  Vec3(float x, float y, float z);

  Vec3();

  Vec3& operator+=(const Vec3 other);

  Vec3 cross(const Vec3 other) const;

  float getMagnitude() const;

  float getMagnitudeSquared() const;

  char* toString(char* singleBuf,
                 std::size_t singleBufSize,
                 char* vecBuf,
                 std::size_t vecBufSize) const;
};

Vec3 operator+(const Vec3& a, const Vec3& b);

Vec3 operator-(const Vec3& a, const Vec3& b);

Vec3 operator*(float s, const Vec3& a);

Vec3 operator/(const Vec3& a, float s);
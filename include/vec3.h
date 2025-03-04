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

  float getMagnitude() const;

  float getMagnitudeSquared() const;

  std::string toString() const;
};

Vec3 operator+(const Vec3& a, const Vec3& b);

Vec3 operator-(const Vec3& a, const Vec3& b);

Vec3 operator*(float s, const Vec3& a);

Vec3 operator/(const Vec3& a, float s);
#pragma once
#include <cmath>
#include <string>

class Vec3 {
 public:
  double x;
  double y;
  double z;

  Vec3(double x, double y, double z);

  Vec3();

  Vec3& operator+=(const Vec3 other);

  double getMagnitude();

  double getMagnitudeSquared();

  std::string toString();
};

Vec3 operator+(const Vec3& a, const Vec3& b);

Vec3 operator-(const Vec3& a, const Vec3& b);

Vec3 operator*(double s, const Vec3& a);
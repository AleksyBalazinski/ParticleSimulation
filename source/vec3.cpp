#include "vec3.h"

Vec3 operator+(const Vec3& a, const Vec3& b) {
  return Vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

Vec3 operator-(const Vec3& a, const Vec3& b) {
  return Vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

Vec3 operator*(double s, const Vec3& a) {
  return Vec3(s * a.x, s * a.y, s * a.z);
}

Vec3 operator/(const Vec3& a, double s) {
  return Vec3(a.x / s, a.y / s, a.z / s);
}

Vec3::Vec3(double x, double y, double z) : x(x), y(y), z(z) {}

Vec3::Vec3() : x(0), y(0), z(0) {}

Vec3& Vec3::operator+=(const Vec3 other) {
  x += other.x;
  y += other.y;
  z += other.z;

  return *this;
}

double Vec3::getMagnitude() const {
  return std::sqrt(x * x + y * y + z * z);
}

double Vec3::getMagnitudeSquared() const {
  return x * x + y * y + z * z;
}

std::string Vec3::toString() const {
  return std::to_string(x) + ' ' + std::to_string(y) + ' ' + std::to_string(z);
}

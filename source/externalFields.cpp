#include "externalFields.h"
#include <cmath>

Vec3 sphRadDecrField(Vec3 pos, Vec3 center, float R, float M, float G) {
  float r = (pos - center).getMagnitude();
  Vec3 dir = (pos - center) / r;
  float g;
  if (r > R) {
    g = -G * M / (r * r);
  } else {
    g = -(G * M / std::powf(R, 3)) * r * (4 - 3 * r / R);
  }

  return g * dir;
}

float sphRadDecrFieldPotential(Vec3 pos, Vec3 center, float R, float M, float G) {
  float r = (pos - center).getMagnitude();
  if (r > R) {
    return -G * M / r;
  }
  float u = r / R;
  return G * M / R * (-2 + u * u * (2 - u));
}

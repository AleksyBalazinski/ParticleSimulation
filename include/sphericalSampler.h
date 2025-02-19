#pragma once

#include <cmath>
#include <random>
#include "vec3.h"

class SphericalSampler {
 private:
  std::default_random_engine re;

  Vec3 getVelocity(Vec3 pos, Vec3 center, float rb, float mb, float rd, float md, float G);

  float sampleFromLinear(float rd);

  inline float implicitInvCDFLinear(float r, float rd, float cdf) {
    return 3 * std::powf(r, 4) - 4 * rd * std::powf(r, 3) + cdf * std::powf(rd, 4);
  }

  inline float implicitInvCDFLinearDer(float r, float rd) {
    return 12 * std::powf(r, 3) - 12 * rd * std::powf(r, 2);
  }

 public:
  SphericalSampler();

  std::vector<Vec3> sample(Vec3 center, float rb, float mb, float rd, float md, float G, int n);
};

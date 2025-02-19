#pragma once

#include <random>
#include "vec3.h"

class DiskSamplerLinear {
 private:
  std::default_random_engine re;

  Vec3 getVelocity(Vec3 pos, Vec3 center, float rb, float mb, float rd, float md, float G);

  float sampleFromLinear(float rd);

  inline float implicitInvCDFLinear(float r, float rd, float cdf) {
    return 2 * std::pow(r, 3) - 3 * rd * r * r + cdf * std::pow(rd, 3);
  }

  inline float implicitInvCDFLinearDer(float r, float rd) { return 6 * r * r - 6 * rd * r; }

 public:
  DiskSamplerLinear();

  std::vector<Vec3>
  sample(Vec3 center, float rb, float mb, float rd, float md, float thickness, float G, int n);
};
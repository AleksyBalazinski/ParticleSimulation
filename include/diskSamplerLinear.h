#pragma once

#include <random>
#include "vec3.h"

class DiskSamplerLinear {
 public:
  DiskSamplerLinear();
  DiskSamplerLinear(unsigned int seed);

  std::vector<Vec3> sample(Vec3 center,
                           float rb,
                           float mb,
                           float rd,
                           float md,
                           float thickness,
                           float G,
                           int n,
                           float r0 = 0.0f);

 private:
  Vec3 getVelocity(Vec3 pos, Vec3 center, float rb, float mb, float rd, float md, float G);

  float sampleFromLinear(float rd, float r0);

  inline float implicitInvCDFLinear(float r, float rd, float cdf, float r0) {
    return 2 * std::powf(r, 3) - 3 * rd * r * r + cdf * std::powf(rd - r0, 2) * (2 * r0 + rd) -
           2 * std::powf(r0, 3) + 3 * rd * r0 * r0;
  }

  inline float implicitInvCDFLinearDer(float r, float rd) { return 6 * r * r - 6 * rd * r; }

  std::default_random_engine re;
};
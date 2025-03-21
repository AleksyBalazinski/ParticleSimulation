#pragma once

#include <random>
#include "external_fields.cuh"
#include "vec3.cuh"

class DiskSamplerLinear {
 private:
  std::default_random_engine re;

  Vec3 getVelocity(Vec3 pos, SphRadDecrFieldParams bulgeParams, float rd, float md, float G);

  float sampleFromLinear(float rd);

  inline float implicitInvCDFLinear(float r, float rd, float cdf) {
    return 2 * std::pow(r, 3) - 3 * rd * r * r + cdf * std::pow(rd, 3);
  }

  inline float implicitInvCDFLinearDer(float r, float rd) { return 6 * r * r - 6 * rd * r; }

 public:
  DiskSamplerLinear();

  std::vector<Vec3> sample(SphRadDecrFieldParams bulgeParams,
                           float rd,
                           float md,
                           float thickness,
                           float G,
                           int n);
};
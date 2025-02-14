#pragma once

#include <cmath>
#include <random>
#include "vec3.h"

class SphericalSampler {
 private:
  std::default_random_engine re;

  Vec3 getVelocity(Vec3 pos, Vec3 center, double rb, double mb, double rd, double md, double G);

  double sampleFromLinear(double rd);

  inline double implicitInvCDFLinear(double r, double rd, double cdf) {
    return 3 * std::pow(r, 4) - 4 * rd * std::pow(r, 3) + cdf * std::pow(rd, 4);
  }

  inline double implicitInvCDFLinearDer(double r, double rd) {
    return 12 * std::pow(r, 3) - 12 * rd * std::pow(r, 2);
  }

 public:
  SphericalSampler();

  std::vector<Vec3>
  sample(Vec3 center, double rb, double mb, double rd, double md, double G, int n);
};

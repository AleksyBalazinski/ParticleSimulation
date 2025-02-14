#pragma once

#include <random>
#include "vec3.h"

class DiskSamplerLinear {
 private:
  std::default_random_engine re;

  Vec3 getVelocity(Vec3 pos, Vec3 center, double rb, double mb, double rd, double md, double G);

  double sampleFromLinear(double rd);

  inline double implicitInvCDFLinear(double r, double rd, double cdf) {
    return 2 * std::pow(r, 3) - 3 * rd * r * r + cdf * std::pow(rd, 3);
  }

  inline double implicitInvCDFLinearDer(double r, double rd) { return 6 * r * r - 6 * rd * r; }

 public:
  DiskSamplerLinear();

  std::vector<Vec3> sample(Vec3 center,
                           double rb,
                           double mb,
                           double rd,
                           double md,
                           double thickness,
                           double G,
                           int n);
};
#pragma once

#include <random>
#include "vec3.h"

class DiskSampler {
 private:
  std::default_random_engine re;

  Vec3 getVelocity(Vec3 pos, Vec3 center, double rb, double mb, double rd, double md, double G);

 public:
  DiskSampler();

  std::vector<Vec3> sample(Vec3 center,
                           double rb,
                           double mb,
                           double rd,
                           double md,
                           double thickness,
                           double G,
                           int n);
};

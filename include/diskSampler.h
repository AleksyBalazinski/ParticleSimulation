#pragma once

#include <random>
#include "vec3.h"

class DiskSampler {
 public:
  DiskSampler(unsigned int seed);

  std::vector<Vec3>
  sample(Vec3 center, float rb, float mb, float rd, float md, float thickness, float G, int n);

 private:
  Vec3 getVelocity(Vec3 pos, Vec3 center, float rb, float mb, float rd, float md, float G);

  std::default_random_engine re;
};

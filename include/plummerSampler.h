#include <numbers>
#include <random>
#include <vector>
#include "vec3.h"

class PlummerSampler {
 public:
  PlummerSampler();
  PlummerSampler(unsigned int seed);

  std::vector<Vec3> sample(Vec3 center, float a, float M, float G, float rMax, int n);

 private:
  float velocityCDF(float x);
  float velocityPDF(float x);

  Vec3 getPosition(Vec3 center,
                   float a,
                   float rMax,
                   std::uniform_real_distribution<float>& u,
                   float& r);
  Vec3 getVelocity(float r, float a, float M, float G, std::uniform_real_distribution<float>& u);

  std::default_random_engine re;
  const float norm;
};
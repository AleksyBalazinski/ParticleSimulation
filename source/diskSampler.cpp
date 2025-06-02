#include "diskSampler.h"
#include <cmath>
#include <numbers>
#include "externalFields.h"
#include "utils.h"

DiskSampler::DiskSampler(unsigned int seed) : re(seed) {}

std::vector<Vec3> DiskSampler::sample(Vec3 center,
                                      float rb,
                                      float mb,
                                      float rd,
                                      float md,
                                      float thickness,
                                      float G,
                                      int n) {
  std::uniform_real_distribution<float> u(0, 1);
  std::vector<Vec3> state(2 * n);

  for (int i = 0; i < n; ++i) {
    float phi = u(re) * 2 * std::numbers::pi_v<float>;
    float r = 0.97f * rd * std::sqrtf(u(re));

    float x = r * std::cos(phi);
    float y = r * std::sin(phi);
    float z = (thickness / 2) * (2 * u(re) - 1);

    state[i] = center + Vec3(x, y, z);
    state[n + i] = getVelocity(state[i], center, rb, mb, rd, md, G);
  }

  return state;
}

Vec3 DiskSampler::getVelocity(Vec3 pos,
                              Vec3 center,
                              float rb,
                              float mb,
                              float rd,
                              float md,
                              float G) {
  Vec3 rVec = pos - center;
  float r = rVec.getMagnitude();

  Vec3 rVecxy = rVec;
  rVecxy.z = 0;
  float rho = rVecxy.getMagnitude();

  float ra = r / rd;
  float n, kSquared;
  n = kSquared = 4 * rd * r / ((rd + r) * (rd + r));
  float k = std::sqrt(kSquared);
  float K = std::comp_ellint_1(k);
  float E = std::comp_ellint_2(k);
  float Pi = std::comp_ellint_3(n, k);

  float density = md / (std::numbers::pi_v<float> * rd * rd);
  float term1 = -(1 + ra * ra) / (ra * (1 + ra)) * K;
  float term2 = (1 + ra) * (2 + ra) / (2 * ra) * E;
  float term3 = -(1 - ra) * (1 - ra) / (2 * (1 + ra)) * Pi;
  float gdVal = 2 * G * density * (term1 + term2 + term3);

  Vec3 gb = sphRadDecrField(pos, center, rb, mb, G);
  Vec3 gd = gdVal * rVecxy / rVecxy.getMagnitude();
  Vec3 g = gb + gd;
  float gVal = g.getMagnitude();

  float v = std::sqrt(gVal * rho * rho / r);
  return Vec3(-v * rVec.y / rho, v * rVec.x / rho, 0);
}

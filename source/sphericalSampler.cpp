#include "sphericalSampler.h"
#include <numbers>
#include "utils.h"

SphericalSampler::SphericalSampler() : re(std::random_device{}()) {}

std::vector<Vec3>
SphericalSampler::sample(Vec3 center, float rb, float mb, float rd, float md, float G, int n) {
  std::uniform_real_distribution<float> u(0, 1);
  std::vector<Vec3> state(2 * n);

  for (int i = 0; i < n; ++i) {
    float phi = u(re) * 2 * std::numbers::pi_v<float>;
    float costheta = u(re) * 2 - 1;
    float theta = std::acos(costheta);
    float r = sampleFromLinear(rd);

    float x = r * std::sin(theta) * std::cos(phi);
    float y = r * std::sin(theta) * std::sin(phi);
    float z = r * std::cos(theta);

    state[i] = center + Vec3(x, y, z);
    state[n + i] = getVelocity(state[i], center, rb, mb, rd, md, G);
  }

  return state;
}

Vec3 SphericalSampler::getVelocity(Vec3 pos,
                                   Vec3 center,
                                   float rb,
                                   float mb,
                                   float rd,
                                   float md,
                                   float G) {
  Vec3 gb = externalFieldBulge(pos, center, rb, mb, G);
  Vec3 gd = externalFieldBulge(pos, center, rd, md, G);
  Vec3 g = gb + gd;
  float gVal = g.getMagnitude();

  Vec3 rVec = pos - center;
  float r = rVec.getMagnitude();

  Vec3 rVecxy = rVec;
  rVecxy.z = 0;
  float rho = rVecxy.getMagnitude();

  float v = std::sqrt(gVal * rho * rho / r);
  return Vec3(-v * rVec.y / rho, v * rVec.x / rho, 0);
}

float SphericalSampler::sampleFromLinear(float rd) {
  std::uniform_real_distribution<float> u(0, 1);

  float cdf = u(re);
  auto f = [rd, cdf, this](float r) { return implicitInvCDFLinear(r, rd, cdf); };
  auto df = [rd, this](float r) { return implicitInvCDFLinearDer(r, rd); };
  return newton(f, df, rd / 2.0f, 100, 0.01f);
}

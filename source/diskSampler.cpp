#include "diskSampler.h"
#include <cmath>
#include <numbers>
#include "utils.h"

Vec3 DiskSampler::getVelocity(Vec3 pos,
                              Vec3 center,
                              double rb,
                              double mb,
                              double rd,
                              double md,
                              double G) {
  Vec3 rVec = pos - center;
  double r = rVec.getMagnitude();

  Vec3 rVecxy = rVec;
  rVecxy.z = 0;
  double rho = rVecxy.getMagnitude();

  double ra = r / rd;
  double n, kSquared;
  n = kSquared = 4 * rd * r / ((rd + r) * (rd + r));
  double k = std::sqrt(kSquared);
  double K = std::comp_ellint_1(k);
  double E = std::comp_ellint_2(k);
  double Pi = std::comp_ellint_3(n, k);

  double density = md / (std::numbers::pi * rd * rd);
  double term1 = -(1 + ra * ra) / (ra * (1 + ra)) * K;
  double term2 = (1 + ra) * (2 + ra) / (2 * ra) * E;
  double term3 = -(1 - ra) * (1 - ra) / (2 * (1 + ra)) * Pi;
  double gdVal = 2 * G * density * (term1 + term2 + term3);

  Vec3 gb = externalFieldBulge(pos, center, rb, mb, G);
  Vec3 gd = gdVal * rVecxy / rVecxy.getMagnitude();
  Vec3 g = gb + gd;
  double gVal = g.getMagnitude();

  double v = std::sqrt(gVal * rho * rho / r);
  return Vec3(-v * rVec.y / rho, v * rVec.x / rho, 0);
}

DiskSampler::DiskSampler() : re(std::random_device{}()) {}

std::vector<Vec3> DiskSampler::sample(Vec3 center,
                                      double rb,
                                      double mb,
                                      double rd,
                                      double md,
                                      double thickness,
                                      double G,
                                      int n) {
  std::uniform_real_distribution<double> u(0, 1);
  std::vector<Vec3> state(2 * n);

  for (int i = 0; i < n; ++i) {
    double phi = u(re) * 2 * std::numbers::pi;
    double r = 0.97 * rd * std::sqrt(u(re));

    double x = r * std::cos(phi);
    double y = r * std::sin(phi);
    double z = (thickness / 2) * (2 * u(re) - 1);

    state[i] = center + Vec3(x, y, z);
    state[n + i] = getVelocity(state[i], center, rb, mb, rd, md, G);
  }

  return state;
}

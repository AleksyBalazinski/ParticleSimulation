#include "diskSamplerLinear.h"
#include <numbers>
#include "utils.h"

Vec3 DiskSamplerLinear::getVelocity(Vec3 pos,
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
  double sigma0 = 3 * md / (std::numbers::pi * rd * rd);
  double k = 2.5;
  double h = 0.66;
  double a = -k / (h * h);
  double gdVal = -G * sigma0 * (a * (ra - h) * (ra - h) + k);

  Vec3 gb = externalFieldBulge(pos, center, rb, mb, G);
  Vec3 gd = gdVal * rVecxy / rVecxy.getMagnitude();
  Vec3 g = gb + gd;
  double gVal = g.getMagnitude();

  double v = std::sqrt(gVal * rho * rho / r);
  return Vec3(-v * rVec.y / rho, v * rVec.x / rho, 0);
}

double DiskSamplerLinear::sampleFromLinear(double rd) {
  std::uniform_real_distribution<double> u(0, 1);

  double cdf = u(re);
  auto f = [rd, cdf, this](double r) { return implicitInvCDFLinear(r, rd, cdf); };
  auto df = [rd, this](double r) { return implicitInvCDFLinearDer(r, rd); };
  return newton(f, df, rd / 2.0, 100, 0.01);
}

DiskSamplerLinear::DiskSamplerLinear() : re(std::random_device{}()) {}

std::vector<Vec3> DiskSamplerLinear::sample(Vec3 center,
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
    double r = sampleFromLinear(rd);

    double x = r * std::cos(phi);
    double y = r * std::sin(phi);
    double z = (thickness / 2) * (2 * u(re) - 1);

    state[i] = center + Vec3(x, y, z);
    state[n + i] = getVelocity(state[i], center, rb, mb, rd, md, G);
  }

  return state;
}

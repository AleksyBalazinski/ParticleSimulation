#include "sphericalSampler.h"
#include <numbers>
#include "utils.h"

double SphericalSampler::newton(std::function<double(double)> f,
                                std::function<double(double)> df,
                                double guess,
                                int maxIter,
                                double tolerance) {
  int i = 0;
  double err = tolerance + 1;
  double x = guess;
  double xPrev;
  while (err > tolerance && i < maxIter) {
    ++i;
    xPrev = x;
    x = x - f(x) / df(x);
    err = std::abs(x - xPrev);
  }

  return x;
}

SphericalSampler::SphericalSampler() : re(std::random_device{}()) {}

std::vector<Vec3>
SphericalSampler::sample(Vec3 center, double rb, double mb, double rd, double md, double G, int n) {
  std::uniform_real_distribution<double> u(0, 1);
  std::vector<Vec3> state(2 * n);

  for (int i = 0; i < n; ++i) {
    double phi = u(re) * 2 * std::numbers::pi;
    double costheta = u(re) * 2 - 1;
    double theta = std::acos(costheta);
    double r = sampleFromLinear(rd);

    double x = r * std::sin(theta) * std::cos(phi);
    double y = r * std::sin(theta) * std::sin(phi);
    double z = r * std::cos(theta);

    state[i] = center + Vec3(x, y, z);
    state[n + i] = getVelocity(state[i], center, rb, mb, rd, md, G);
  }

  return state;
}

Vec3 SphericalSampler::getVelocity(Vec3 pos,
                                   Vec3 center,
                                   double rb,
                                   double mb,
                                   double rd,
                                   double md,
                                   double G) {
  Vec3 gb = externalFieldBulge(pos, center, rb, mb, G);
  Vec3 gd = externalFieldBulge(pos, center, rd, md, G);
  Vec3 g = gb + gd;
  double gVal = g.getMagnitude();

  Vec3 rVec = pos - center;
  double r = rVec.getMagnitude();

  Vec3 rVecxy = rVec;
  rVecxy.z = 0;
  double rho = rVecxy.getMagnitude();

  double v = std::sqrt(gVal * rho * rho / r);
  return Vec3(-v * rVec.y / rho, v * rVec.x / rho, 0);
}

double SphericalSampler::sampleFromLinear(double rd) {
  std::uniform_real_distribution<double> u(0, 1);

  double cdf = u(re);
  auto f = [rd, cdf, this](double r) { return implicitInvCDFLinear(r, rd, cdf); };
  auto df = [rd, this](double r) { return implicitInvCDFLinearDer(r, rd); };
  return newton(f, df, rd / 2.0, 100, 0.01);
}

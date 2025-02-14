#include "utils.h"
#include <random>

std::vector<Vec3> randomInitialState(int particlesCnt,
                                     std::pair<Vec3, Vec3> initalPositionRange,
                                     std::pair<Vec3, Vec3> initialVelocityRange) {
  std::default_random_engine re;

  std::uniform_real_distribution<double> posDistX(initalPositionRange.first.x,
                                                  initalPositionRange.second.x);
  std::uniform_real_distribution<double> posDistY(initalPositionRange.first.y,
                                                  initalPositionRange.second.y);
  std::uniform_real_distribution<double> posDistZ(initalPositionRange.first.z,
                                                  initalPositionRange.second.z);

  std::uniform_real_distribution<double> velocityDistX(initialVelocityRange.first.x,
                                                       initialVelocityRange.second.x);
  std::uniform_real_distribution<double> velocityDistY(initialVelocityRange.first.y,
                                                       initialVelocityRange.second.y);
  std::uniform_real_distribution<double> velocityDistZ(initialVelocityRange.first.z,
                                                       initialVelocityRange.second.z);

  std::vector<Vec3> initialState(2 * particlesCnt);

  for (int i = 0; i < particlesCnt; i++) {
    initialState[i].x = posDistX(re);
    initialState[i].y = posDistY(re);
    initialState[i].z = posDistZ(re);

    initialState[particlesCnt + i].x = velocityDistX(re);
    initialState[particlesCnt + i].y = velocityDistY(re);
    initialState[particlesCnt + i].z = velocityDistZ(re);
  }

  return initialState;
}

std::vector<double> randomMasses(int particlesCnt, std::pair<double, double> massRange) {
  std::default_random_engine re;
  std::uniform_real_distribution<double> massDist(massRange.first, massRange.second);

  std::vector<double> masses(particlesCnt);

  for (int i = 0; i < particlesCnt; i++) {
    masses[i] = massDist(re);
  }

  return masses;
}

Vec3 externalFieldBulge(Vec3 pos, Vec3 bulge, double rb, double mb, double G) {
  double r = (pos - bulge).getMagnitude();
  Vec3 dir = (pos - bulge) / r;
  double g;
  if (r > rb) {
    g = -G * mb / (r * r);
  } else {
    g = -(G * mb / std::pow(rb, 3)) * r * (4 - 3 * r / rb);
  }

  return g * dir;
}

double newton(std::function<double(double)> f,
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
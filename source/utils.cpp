#include "utils.h"
#include <random>

std::vector<Vec3> randomInitialState(int particlesCnt,
                                     std::pair<Vec3, Vec3> initalPositionRange,
                                     std::pair<Vec3, Vec3> initialVelocityRange) {
  std::default_random_engine re;

  std::uniform_real_distribution<float> posDistX(initalPositionRange.first.x,
                                                 initalPositionRange.second.x);
  std::uniform_real_distribution<float> posDistY(initalPositionRange.first.y,
                                                 initalPositionRange.second.y);
  std::uniform_real_distribution<float> posDistZ(initalPositionRange.first.z,
                                                 initalPositionRange.second.z);

  std::uniform_real_distribution<float> velocityDistX(initialVelocityRange.first.x,
                                                      initialVelocityRange.second.x);
  std::uniform_real_distribution<float> velocityDistY(initialVelocityRange.first.y,
                                                      initialVelocityRange.second.y);
  std::uniform_real_distribution<float> velocityDistZ(initialVelocityRange.first.z,
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

std::vector<float> randomMasses(int particlesCnt, std::pair<float, float> massRange) {
  std::default_random_engine re;
  std::uniform_real_distribution<float> massDist(massRange.first, massRange.second);

  std::vector<float> masses(particlesCnt);

  for (int i = 0; i < particlesCnt; i++) {
    masses[i] = massDist(re);
  }

  return masses;
}

Vec3 externalFieldBulge(Vec3 pos, Vec3 bulge, float rb, float mb, float G) {
  float r = (pos - bulge).getMagnitude();
  Vec3 dir = (pos - bulge) / r;
  float g;
  if (r > rb) {
    g = -G * mb / (r * r);
  } else {
    g = -(G * mb / std::powf(rb, 3)) * r * (4 - 3 * r / rb);
  }

  return g * dir;
}

float newton(std::function<float(float)> f,
             std::function<float(float)> df,
             float guess,
             int maxIter,
             float tolerance) {
  int i = 0;
  float err = tolerance + 1;
  float x = guess;
  float xPrev;
  while (err > tolerance && i < maxIter) {
    ++i;
    xPrev = x;
    x = x - f(x) / df(x);
    err = std::abs(x - xPrev);
  }

  return x;
}
#include "plummerSampler.h"
#include <iostream>
#include <numbers>
#include "utils.h"

PlummerSampler::PlummerSampler() : PlummerSampler(std::random_device{}()) {}

PlummerSampler::PlummerSampler(unsigned int seed)
    : re(seed), norm(512.0f / (7 * std::numbers::pi_v<float>)) {}

std::vector<Vec3> PlummerSampler::sample(Vec3 center,
                                         float a,
                                         float rMax,
                                         float M,
                                         float G,
                                         int n) {
  std::vector<Vec3> state(2 * n);
  std::uniform_real_distribution<float> u(0, 1);

  for (int i = 0; i < n; ++i) {
    float r;
    state[i] = getPosition(center, a, rMax, u, r);
    state[n + i] = getVelocity(r, a, M, G, u);
  }

  return state;
}

float PlummerSampler::velocityCDF(float x) {
  float sqrtTerm = std::sqrtf(1 - x * x);
  float arctanTerm = std::atanf(x / sqrtTerm);
  float numerator = (x * sqrtTerm *
                         (-105 + 1210 * x * x - 2104 * std::powf(x, 4) + 1488 * std::powf(x, 6) -
                          384 * std::powf(x, 8)) +
                     105 * arctanTerm);
  return (numerator / 3840.0f) * norm;
}

float PlummerSampler::velocityPDF(float x) {
  return x * x * std::powf(1 - x * x, 3.5f) * norm;
}

Vec3 PlummerSampler::getPosition(Vec3 center,
                                 float a,
                                 float rMax,
                                 std::uniform_real_distribution<float>& u,
                                 float& rOut) {
  float phi = 2 * std::numbers::pi_v<float> * u(re);
  float r = a * std::powf(std::powf(u(re), -2.0f / 3) - 1, -1.0 / 2);
  if (r > rMax) {
    r = rMax;
  }
  float theta = std::acosf(1 - 2 * u(re));

  float x = r * std::sinf(theta) * std::cosf(phi);
  float y = r * std::sinf(theta) * std::sinf(phi);
  float z = r * std::cosf(theta);

  rOut = r;
  return center + Vec3(x, y, z);
}

Vec3 PlummerSampler::getVelocity(float r,
                                 float a,
                                 float M,
                                 float G,
                                 std::uniform_real_distribution<float>& u) {
  float potential = -G * M / std::sqrtf(r * r + a * a);
  float vEsc = std::sqrtf(-2 * potential);
  float cdf = u(re);
  auto f = [this, cdf](float x) -> float { return velocityCDF(x) - cdf; };
  auto df = std::bind(&PlummerSampler::velocityPDF, this, std::placeholders::_1);
  float q = newton(f, df, 0.5f, 100, 1e-3f);
  float v = q * vEsc;

  float phi = 2 * std::numbers::pi_v<float> * u(re);
  float theta = std::acosf(1 - 2 * u(re));
  float vx = v * std::sinf(theta) * std::cosf(phi);
  float vy = v * std::sinf(theta) * std::sinf(phi);
  float vz = v * std::cosf(theta);

  return Vec3(vx, vy, vz);
}

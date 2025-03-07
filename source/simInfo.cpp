#include "simInfo.h"

float totalEnergy(const std::vector<Vec3>& state, const std::vector<float>& masses, float G) {
  const int n = (int)masses.size();
  float potentialEnergy = 0;
  float kineticEnergy = 0;
  for (int i = 0; i < n; i++) {
    kineticEnergy += 0.5f * masses[i] * state[n + i].getMagnitudeSquared();
    for (int j = i + 1; j < n; j++) {
      auto rij = state[i] - state[j];
      potentialEnergy += (-1) * G * masses[i] * masses[j] / rij.getMagnitude();
    }
  }

  return potentialEnergy + kineticEnergy;
}

float potentialEnergy(std::vector<Vec3>::iterator posBegin,
                      std::vector<Vec3>::iterator posEnd,
                      const std::vector<float>& masses,
                      float G) {
  const int n = (int)masses.size();
  float potentialEnergy = 0;
  auto posIt = posBegin;
  for (int i = 0; posIt != posEnd; ++posIt, ++i) {
    auto posIt2 = posIt + 1;
    for (int j = i + 1; posIt2 != posEnd; ++posIt2, ++j) {
      auto rij = *posIt - *posIt2;
      potentialEnergy += (-1) * G * masses[i] * masses[j] / rij.getMagnitude();
    }
  }

  return potentialEnergy;
}

float kineticEnergy(std::vector<Vec3>::iterator vBegin,
                    std::vector<Vec3>::iterator vEnd,
                    const std::vector<float>& masses,
                    float G) {
  const int n = (int)masses.size();
  float kineticEnergy = 0;

  auto vIt = vBegin;
  for (int i = 0; vIt != vEnd; ++vIt, ++i) {
    kineticEnergy += 0.5f * masses[i] * vIt->getMagnitudeSquared();
  }

  return kineticEnergy;
}

float totalEnergy(std::vector<Vec3>::iterator posBegin,
                  std::vector<Vec3>::iterator posEnd,
                  std::vector<Vec3>::iterator vBegin,
                  std::vector<Vec3>::iterator vEnd,
                  const std::vector<float>& masses,
                  float G) {
  const int n = (int)masses.size();
  float potentialEnergy = 0;
  float kineticEnergy = 0;

  auto posIt = posBegin;
  auto vIt = vBegin;
  for (int i = 0; posIt != posEnd && vIt != vEnd; ++posIt, ++vIt, ++i) {
    kineticEnergy += 0.5f * masses[i] * vIt->getMagnitudeSquared();
    auto posIt2 = posIt + 1;
    for (int j = i + 1; posIt2 != posEnd; ++posIt2, ++j) {
      auto rij = *posIt - *posIt2;
      potentialEnergy += (-1) * G * masses[i] * masses[j] / rij.getMagnitude();
    }
  }

  return potentialEnergy + kineticEnergy;
}

Vec3 totalMomentum(const std::vector<Vec3>& state, const std::vector<float>& masses) {
  const int n = (int)masses.size();
  Vec3 momentum;
  for (int i = 0; i < n; i++) {
    momentum += masses[i] * state[n + i];
  }
  return momentum;
}

Vec3 totalMomentum(std::vector<Vec3>::iterator vBegin,
                   std::vector<Vec3>::iterator vEnd,
                   const std::vector<float>& masses) {
  Vec3 momentum;
  auto vIt = vBegin;
  for (int i = 0; vIt != vEnd; ++vIt, ++i) {
    momentum += masses[i] * (*vIt);
  }
  return momentum;
}

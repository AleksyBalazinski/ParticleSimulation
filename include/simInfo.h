#pragma once

#include <vector>
#include "particle.h"
#include "vec3.h"

class SimInfo {
 private:
  Vec3 expectedMomentum;

 public:
  static float potentialEnergy(std::vector<Vec3>::iterator posBegin,
                               std::vector<Vec3>::iterator posEnd,
                               const std::vector<float>& masses,
                               float G);

  static float potentialEnergy(const std::vector<Particle>& particles, float G);

  static float kineticEnergy(std::vector<Vec3>::iterator vBegin,
                             std::vector<Vec3>::iterator vEnd,
                             const std::vector<float>& masses,
                             float G);

  static float kineticEnergy(const std::vector<Particle>& particles, float G);

  static float totalEnergy(const std::vector<Vec3>& state,
                           const std::vector<float>& masses,
                           float G);

  static float totalEnergy(std::vector<Vec3>::iterator posBegin,
                           std::vector<Vec3>::iterator posEnd,
                           std::vector<Vec3>::iterator vBegin,
                           std::vector<Vec3>::iterator vEnd,
                           const std::vector<float>& masses,
                           float G);

  static Vec3 totalMomentum(const std::vector<Vec3>& state, const std::vector<float>& masses);

  static Vec3 totalMomentum(std::vector<Vec3>::iterator vBegin,
                            std::vector<Vec3>::iterator vEnd,
                            const std::vector<float>& masses);

  static Vec3 totalMomentum(const std::vector<Particle>& particles);

  void setInitialMomentum(const std::vector<Particle>& particles);

  Vec3 updateExpectedMomentum(Vec3 externalForce, float DT);
};

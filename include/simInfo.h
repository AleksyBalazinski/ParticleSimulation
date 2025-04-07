#pragma once

#include <functional>
#include <vector>
#include "grid.h"
#include "particle.h"
#include "vec3.h"

class SimInfo {
 public:
  static float potentialEnergy(std::vector<Vec3>::iterator posBegin,
                               std::vector<Vec3>::iterator posEnd,
                               const std::vector<float>& masses,
                               float G);
  static float kineticEnergy(std::vector<Vec3>::iterator vBegin,
                             std::vector<Vec3>::iterator vEnd,
                             const std::vector<float>& masses,
                             float G);

  static Vec3 totalMomentum(std::vector<Vec3>::iterator vBegin,
                            std::vector<Vec3>::iterator vEnd,
                            const std::vector<float>& masses);

  static float potentialEnergy(const Grid& grid,
                               const std::vector<Particle>& particles,
                               std::function<float(Vec3)> externalPotential,
                               float H,
                               float DT,
                               float G);
  static float potentialEnergy(const std::vector<std::complex<float>>& gridDensity,
                               const std::vector<std::complex<float>>& gridPotential,
                               const std::vector<Particle>& particles,
                               std::function<float(Vec3)> externalPotential,
                               float H,
                               float DT,
                               float G);

  static float kineticEnergy(const std::vector<Particle>& particles);

  static Vec3 totalMomentum(const std::vector<Particle>& particles);
  static Vec3 totalAngularMomentum(const std::vector<Particle>& particles);

  void setInitialMomentum(const std::vector<Particle>& particles);
  Vec3 updateExpectedMomentum(Vec3 externalForce, float DT);

 private:
  Vec3 expectedMomentum;
};

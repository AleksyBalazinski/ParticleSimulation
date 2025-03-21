#include "simInfo.h"
#include "unitConversions.h"

float SimInfo::potentialEnergy(std::vector<Vec3>::iterator posBegin,
                               std::vector<Vec3>::iterator posEnd,
                               const std::vector<float>& masses,
                               float G) {
  const int n = (int)masses.size();
  float potentialEnergy = 0;
  auto posIt = posBegin;
  float eps = 0.01f;
  for (int i = 0; posIt != posEnd; ++posIt, ++i) {
    auto posIt2 = posIt + 1;
    for (int j = i + 1; posIt2 != posEnd; ++posIt2, ++j) {
      auto rij = *posIt - *posIt2;
      potentialEnergy +=
          (-1) * G * masses[i] * masses[j] / std::powf(rij.getMagnitudeSquared() + eps * eps, 0.5f);
    }
  }

  return potentialEnergy;
}

float SimInfo::kineticEnergy(std::vector<Vec3>::iterator vBegin,
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

Vec3 SimInfo::totalMomentum(std::vector<Vec3>::iterator vBegin,
                            std::vector<Vec3>::iterator vEnd,
                            const std::vector<float>& masses) {
  Vec3 momentum;
  auto vIt = vBegin;
  for (int i = 0; vIt != vEnd; ++vIt, ++i) {
    momentum += masses[i] * (*vIt);
  }
  return momentum;
}

float SimInfo::potentialEnergy(const Grid& grid,
                               const std::vector<Particle>& particles,
                               std::function<float(Vec3)> externalPotential,
                               float H,
                               float DT,
                               float G) {
  float internal = 0;
  for (int i = 0; i < grid.getLength(); ++i) {
    auto [x, y, z] = grid.indexTripleFromFlat(i);
    internal += densityToOriginalUnits(grid.getDensity(x, y, z), DT, G) *
                potentialToOriginalUnits(grid.getPotential(x, y, z), H, DT);
  }

  float external = 0;
  for (const auto& p : particles) {
    external += p.mass * externalPotential(p.position);
  }

  float vol = H * H * H;
  return 0.5f * vol * internal + external;
}

float SimInfo::kineticEnergy(const std::vector<Particle>& particles) {
  float ke = 0;
  for (const auto& p : particles) {
    ke += 0.5f * p.mass * p.velocity.getMagnitudeSquared();
  }

  return ke;
}

Vec3 SimInfo::totalMomentum(const std::vector<Particle>& particles) {
  Vec3 momentum;
  for (const auto& p : particles) {
    momentum += p.mass * p.integerStepVelocity;
  }
  return momentum;
}

Vec3 SimInfo::totalAngularMomentum(const std::vector<Particle>& particles) {
  Vec3 angularMomentum;
  for (const auto& p : particles) {
    angularMomentum += p.mass * p.position.cross(p.integerStepVelocity);
  }

  return angularMomentum;
}

void SimInfo::setInitialMomentum(const std::vector<Particle>& particles) {
  expectedMomentum = totalMomentum(particles);
}

Vec3 SimInfo::updateExpectedMomentum(Vec3 externalForce, float DT) {
  expectedMomentum += DT * externalForce;
  return expectedMomentum;
}

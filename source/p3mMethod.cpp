#include "p3mMethod.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <numbers>
#include "leapfrog.h"
#include "unitConversions.h"

P3MMethod::P3MMethod(PMMethod& pmMethod,
                     float compBoxSize,
                     float cutoffRadius,
                     float particleDiameter,
                     float H,
                     bool useSRForceTable)
    : pmMethod(pmMethod),
      chainingMesh(compBoxSize, cutoffRadius, H),
      cutoffRadius(lengthToCodeUnits(cutoffRadius, pmMethod.getH())),
      particleDiameter(lengthToCodeUnits(particleDiameter, H)),
      useSRForceTable(useSRForceTable) {
  this->tabulatedValuesCnt = 500;
  float re = this->cutoffRadius;
  this->deltaSquared = re * re / (this->tabulatedValuesCnt - 1);
  if (useSRForceTable) {
    initSRForceTable();
  }
}

void P3MMethod::calculateShortRangeForces(std::vector<Particle>& particles) {
  for (auto& p : particles) {
    p.shortRangeForce = Vec3();
  }

  chainingMesh.clear();
  chainingMesh.fill(particles);

  for (int q = 0; q < chainingMesh.getSize(); ++q) {
    auto neighbors = chainingMesh.getNeighborsAndSelf(q);
    for (int qn : neighbors) {
      if (qn == -1) {
        continue;
      }
      for (auto node = chainingMesh.getParticlesInCell(q); node != nullptr; node = node->next) {
        for (auto nodeN = chainingMesh.getParticlesInCell(qn); nodeN != nullptr;
             nodeN = nodeN->next) {
          updateShortRangeFoces(node->particleId, nodeN->particleId, q, qn, particles);
        }
      }
    }
  }
}

void correctAccelerations(std::vector<Particle>& particles) {
  for (auto& p : particles) {
    p.acceleration += p.shortRangeForce / p.mass;
  }
}

void P3MMethod::run(const int simLength,
                    bool collectDiagnostics,
                    bool recordField,
                    const char* positionsPath,
                    const char* energyPath,
                    const char* momentumPath,
                    const char* expectedMomentumPath,
                    const char* fieldPath) {
  std::vector<Particle>& particles = pmMethod.getParticles();
  StateRecorder stateRecorder(positionsPath, energyPath, momentumPath, expectedMomentumPath,
                              fieldPath);
  float H = pmMethod.getH();
  float DT = pmMethod.getDT();
  float G = pmMethod.getG();

  if (collectDiagnostics) {
    simInfo.setInitialMomentum(particles);
  }

  stateToCodeUnits(particles, H, DT);
  massToCodeUnits(particles, H, DT, G);

  pmMethod.initGreensFunction();
  pmMethod.pmMethodStep();
  calculateShortRangeForces(particles);
  correctAccelerations(particles);

  setHalfStepVelocities(particles);

  for (int t = 0; t <= simLength; ++t) {
    std::cout << "progress: " << float(t) / simLength << '\r';
    std::cout.flush();

    updatePositions(particles);

    if (collectDiagnostics) {
      setIntegerStepVelocities(particles);
      integerStepVelocitiesToOriginalUnits(particles, H, DT);
    }

    stateToOriginalUnits(particles, H, DT);

    stateRecorder.recordPositions(particles);

    if (collectDiagnostics) {
      massToOriginalUnits(particles, H, DT, G);
      auto expectedMomentum =
          simInfo.updateExpectedMomentum(pmMethod.totalExternalForceOrigUnits(), DT);
      stateRecorder.recordExpectedMomentum(expectedMomentum);
      stateRecorder.recordTotalMomentum(SimInfo::totalMomentum(particles));
      auto pe = SimInfo::potentialEnergy(pmMethod.getGrid(), particles,
                                         pmMethod.getExternalPotential(), H, DT, G);
      auto ke = SimInfo::kineticEnergy(particles);
      stateRecorder.recordEnergy(pe, ke);
    }
    if (recordField) {
      stateRecorder.recordField(particles, H, DT);
    }

    if (pmMethod.escapedComputationalBox()) {
      std::cout << "Particle moved outside the computational box.\n";
      break;
    }

    stateToCodeUnits(particles, H, DT);
    if (collectDiagnostics) {
      massToCodeUnits(particles, H, DT, G);
      integerStepVelocitiesToCodeUnits(particles, H, DT);
    }

    pmMethod.pmMethodStep();
    calculateShortRangeForces(particles);
    correctAccelerations(particles);

    updateVelocities(particles);
  }

  stateRecorder.flush();
}

float P3MMethod::referenceForceS1(float r, float a) {
  float G = 1 / (4 * std::numbers::pi_v<float>);
  float eps = 0.01f;
  if (r >= a) {
    return G / (r * r + eps * eps);
  }
  return G / (a * a) * (8 * r / a - 9 * r * r / (a * a) + 2 * std::powf(r / a, 4));
}

Vec3 P3MMethod::shortRangeForce(Vec3 rij, float mi, float mj, float a) {
  float rijLength = rij.getMagnitude();
  Vec3 rijDir = rij / rijLength;
  float G = 1 / (4 * std::numbers::pi_v<float>);

  Vec3 Rij = -mi * mj * referenceForceS1(rijLength, particleDiameter) * rijDir;
  float eps = 0.01f;
  Vec3 totalForceij = -G * mi * mj / (rijLength * rijLength + eps * eps) * rijDir;
  return totalForceij - Rij;
}

Vec3 P3MMethod::shortRangeForceFromTable(Vec3 rij, float mi, float mj, float a) {
  float ksi = rij.getMagnitudeSquared() / deltaSquared;
  int t = int(ksi);
  assert(t < tabulatedValuesCnt - 1);
  float F = mi * mj * (FTable[t] + (ksi - t) * (FTable[t + 1] - FTable[t]));
  return F * rij;
}

void P3MMethod::updateShortRangeFoces(int i,
                                      int j,
                                      int q,
                                      int qn,
                                      std::vector<Particle>& particles) {
  if (i == j) {
    return;
  }

  Vec3 rij = particles[i].position - particles[j].position;
  if (rij.getMagnitudeSquared() >= cutoffRadius * cutoffRadius) {
    return;
  }

  Vec3 shortRangeij;
  if (!useSRForceTable) {
    shortRangeij = shortRangeForce(rij, particles[i].mass, particles[j].mass, particleDiameter);
  } else {
    shortRangeij =
        shortRangeForceFromTable(rij, particles[i].mass, particles[j].mass, particleDiameter);
  }

  particles[i].shortRangeForce += shortRangeij;
  if (qn != q) {
    particles[j].shortRangeForce += -1 * shortRangeij;
  }
}

void P3MMethod::initSRForceTable() {  // TODO: add S2 reference force
  for (int i = 0; i < tabulatedValuesCnt; ++i) {
    float rSquared = i * deltaSquared;
    float r = std::sqrtf(rSquared);
    float R = -referenceForceS1(r, particleDiameter);

    float G = 1 / (4 * std::numbers::pi_v<float>);
    float eps = 0.01f;
    float totalForce = -G / (r * r + eps * eps);

    FTable.push_back(r == 0 ? 0 : totalForce / r);
  }
}

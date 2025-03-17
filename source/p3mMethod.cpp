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
                     float particleRadius,
                     float H)
    : pmMethod(pmMethod),
      chainingMesh(compBoxSize, cutoffRadius, H),
      cutoffRadius(lengthToCodeUnits(cutoffRadius, pmMethod.getH())),
      particleRadius(lengthToCodeUnits(particleRadius, H)) {
  this->tabulatedValuesCnt = 500;
  float re = this->cutoffRadius;
  this->deltaSquared = re * re / (this->tabulatedValuesCnt - 1);
  for (int i = 0; i < this->tabulatedValuesCnt; ++i) {
    float rSquared = i * this->deltaSquared;
    float r = std::sqrtf(rSquared);
    float R = -referenceForce(r, this->particleRadius);

    float G = 1 / (4 * std::numbers::pi_v<float>);
    float eps = 0.05f;
    float totalForce = -G / (r * r + eps * eps);

    FTable.push_back(r == 0 ? 0 : totalForce / r);
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
          int i = node->particleId;
          int j = nodeN->particleId;
          if (i == j) {
            continue;
          }
          Vec3 rij = particles[i].position - particles[j].position;
          if (rij.getMagnitudeSquared() >= cutoffRadius * cutoffRadius) {
            continue;
          }

          float rijLength = rij.getMagnitude();
          Vec3 rijDir = rij / rijLength;
          float mi = particles[i].mass;
          float mj = particles[j].mass;
          float G = 1 / (4 * std::numbers::pi_v<float>);

          Vec3 Rij = -mi * mj * referenceForce(rijLength, particleRadius) * rijDir;
          float eps = 0.05f;
          Vec3 totalForceij = -G * mi * mj / (rijLength * rijLength + eps * eps) * rijDir;
          Vec3 shortRangeij = totalForceij - Rij;

          particles[i].shortRangeForce += shortRangeij;
          if (qn != q) {
            particles[j].shortRangeForce += -1 * shortRangeij;
          }

          // float ksi = rij.getMagnitudeSquared() / deltaSquared;
          // int t = int(ksi);
          // assert(t < tabulatedValuesCnt - 1);
          // float mi = particles[i].mass;
          // float mj = particles[j].mass;
          // float F = mi * mj * (FTable[t] + (ksi - t) * (FTable[t + 1] - FTable[t]));
          // particles[i].shortRangeForce += F * rij;
          // if (qn != q) {
          //   particles[j].shortRangeForce += -F * rij;
          // }
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

float P3MMethod::referenceForce(float r, float a) {
  float G = 1 / (4 * std::numbers::pi_v<float>);
  float eps = 0.05f;
  if (r >= a) {
    return G / (r * r + eps * eps);
  }
  return G / (a * a) * (8 * r / a - 9 * r * r / (a * a) + 2 * std::powf(r / a, 4));
}

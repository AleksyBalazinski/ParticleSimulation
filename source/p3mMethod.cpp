#include "p3mMethod.h"
#include <algorithm>
#include <cmath>
#include <execution>
#include <iostream>
#include <numbers>
#include <numeric>
#include <ranges>
#include <thread>
#include "leapfrog.h"
#include "measureTime.h"
#include "unitConversions.h"

declareTimeAcc(chainingMeshSetup);
declareTimeAcc(shortRangeForcesCalc);
declareTimeAcc(pmStep);
declareTimeAcc(correctAccelerations);

P3MMethod::P3MMethod(
#ifdef CUDA
    PMMethodGPU& pmMethod,
#else
    PMMethod& pmMethod,
#endif
    std::tuple<float, float, float> compBoxSize,
    float cutoffRadius,
    float particleDiameter,
    float H,
    float softeningLength,
    CloudShape cloudShape,
    bool useSRForceTable,
    bool enableYSorting)
    : pmMethod(pmMethod),
      chainingMesh(compBoxSize, cutoffRadius, H, int(pmMethod.getParticles().size())),
      cutoffRadius(lengthToCodeUnits(cutoffRadius, pmMethod.getH())),
      particleDiameter(lengthToCodeUnits(particleDiameter, H)),
      softeningLength(lengthToCodeUnits(softeningLength, H)),
      useSRForceTable(useSRForceTable),
      enableYSorting(enableYSorting),
      cloudShape(cloudShape),
      threadTimes(std::thread::hardware_concurrency(), 0.0f) {
  this->tabulatedValuesCnt = 500;
  float re = this->cutoffRadius;
  this->deltaSquared = re * re / (this->tabulatedValuesCnt - 1);
  if (useSRForceTable) {
    initSRForceTable();
  }
}

void correctAccelerations(std::vector<Particle>& particles) {
  std::for_each(std::execution::par, particles.begin(), particles.end(), [](Particle& p) {
    Vec3 totalSRForce = std::accumulate(p.shortRangeFromNeighbor.begin(),
                                        p.shortRangeFromNeighbor.end(), Vec3::zero()) +
                        p.shortRangeForce;
    p.acceleration += totalSRForce / p.mass;
  });
}

void P3MMethod::run(StateRecorder& stateRecorder,
                    const int simLength,
                    bool collectDiagnostics,
                    bool recordField) {
  std::vector<Particle>& particles = pmMethod.getParticles();

  float H = pmMethod.getH();
  float DT = pmMethod.getDT();
  float G = pmMethod.getG();

  if (collectDiagnostics) {
    simInfo.setInitialMomentum(particles);
  }

  stateToCodeUnits(particles, H, DT);
  massToCodeUnits(particles, H, DT, G);

  pmMethod.initGreensFunction();
#ifdef CUDA
  pmMethod.copyParticlesHostToDevice();
#endif
  measureTime(pmStep, pmMethod.pmMethodStep());
#ifdef CUDA
  pmMethod.copyParticlesDeviceToHost();
#endif
  // auto start = std::chrono::high_resolution_clock::now();
  measureTime(shortRangeForcesCalc, calculateShortRangeForces(particles));
  measureTime(correctAccelerations, correctAccelerations(particles));
  // auto end = std::chrono::high_resolution_clock::now();
  // auto elapsed = duration_cast<std::chrono::milliseconds>(end - start);
  // std::cout << particleDiameter << ", " << elapsed.count() << std::endl;

  setHalfStepVelocities(particles);

  for (int t = 0; t <= simLength; ++t) {
    std::cout << "progress: " << float(t) / simLength << '\r';
    std::cout.flush();

    if (recordField) {
      stateRecorder.recordField(particles, H, DT);
    }

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
      stateRecorder.recordTotalAngularMomentum(SimInfo::totalAngularMomentum(particles));
#ifdef CUDA
      pmMethod.copyGridDensityToHost();
      pmMethod.copyGridPotentialToHost();
      auto pe = SimInfo::potentialEnergy(pmMethod.getGridDensity(), pmMethod.getGridPotential(),
                                         particles, pmMethod.getExternalPotential(), H, DT, G);
#else
      auto pe = SimInfo::potentialEnergy(pmMethod.getGrid(), particles,
                                         pmMethod.getExternalPotential(), H, DT, G);
#endif
      auto ke = SimInfo::kineticEnergy(particles);
      stateRecorder.recordEnergy(pe, ke);
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

#ifdef CUDA
    pmMethod.copyParticlesHostToDevice();
#endif
    measureTime(pmStep, pmMethod.pmMethodStep());
#ifdef CUDA
    pmMethod.copyParticlesDeviceToHost();
#endif
    measureTime(shortRangeForcesCalc, calculateShortRangeForces(particles));
    measureTime(correctAccelerations, correctAccelerations(particles));

    updateVelocities(particles);
  }

  printTime(chainingMeshSetup);
  printTime(shortRangeForcesCalc);
  printTime(pmStep);
  printTime(correctAccelerations);

  for (int i = 0; i < threadTimes.size(); ++i) {
    std::cout << "thread " << i << " was working for " << threadTimes[i] << " ms\n";
  }

  stateRecorder.flush();
}

void P3MMethod::calculateShortRangeForces(std::vector<Particle>& particles) {
  const int maxThreads = std::thread::hardware_concurrency();
  std::for_each(std::execution::par, particles.begin(), particles.end(), [](Particle& p) {
    p.shortRangeForce = Vec3::zero();
    for (auto& srForceNeighbor : p.shortRangeFromNeighbor) {
      srForceNeighbor = Vec3::zero();
    }
  });

  if (enableYSorting) {
    measureTime(chainingMeshSetup, chainingMesh.fillWithYSorting(particles));
  } else {
    measureTime(chainingMeshSetup, chainingMesh.fill(particles));
  }

  std::vector<std::thread> smallCellThreads;
  for (int tid = 0; tid < maxThreads; ++tid) {
    smallCellThreads.emplace_back(&P3MMethod::updateSRForcesThreadJob, this, tid, maxThreads,
                                  std::ref(particles));
  }

  for (auto& t : smallCellThreads) {
    t.join();
  }
}

float P3MMethod::referenceForceS1(float r) {
  const float a = particleDiameter;
  const float G = 1 / (4 * std::numbers::pi_v<float>);
  if (r >= a) {
    return G / (r * r);
  }
  return G / (a * a) * (8 * r / a - 9 * r * r / (a * a) + 2 * std::powf(r / a, 4));
}

float P3MMethod::referenceForceS2(float r) {
  const float a = particleDiameter;
  const float G = 1 / (4 * std::numbers::pi_v<float>);
  const float u = 2 * r / a;
  if (u <= 1) {
    return G / (35 * std::powf(a, 2)) *
           (224 * u - 224 * std::powf(u, 3) + 70 * std::powf(u, 4) + 48 * std::powf(u, 5) -
            21 * std::powf(u, 6));
  }
  if (u <= 2) {
    return G / (35 * std::powf(a, 2)) *
           (12 / std::powf(u, 2) - 224 + 896 * u - 840 * std::powf(u, 2) + 224 * std::powf(u, 3) +
            70 * std::powf(u, 4) - 48 * std::powf(u, 5) + 7 * std::powf(u, 6));
  }
  return G / (r * r);
}

Vec3 P3MMethod::shortRangeForce(Vec3 rij, float mi, float mj) {
  float rijLength = rij.getMagnitude();
  Vec3 rijDir = rij / rijLength;
  float G = 1 / (4 * std::numbers::pi_v<float>);

  float R;
  if (cloudShape == CloudShape::S1) {
    R = referenceForceS1(rijLength);
  } else if (cloudShape == CloudShape::S2) {
    R = referenceForceS2(rijLength);
  } else {
    throw std::invalid_argument("not implemnted");
  }

  Vec3 Rij = -mi * mj * R * rijDir;
  float eps = softeningLength;
  Vec3 totalForceij = -G * mi * mj / (rijLength * rijLength + eps * eps) * rijDir;
  return totalForceij - Rij;
}

Vec3 P3MMethod::shortRangeForceFromTable(Vec3 rij, float mi, float mj) {
  float ksi = rij.getMagnitudeSquared() / deltaSquared;
  int t = int(ksi);
  float F = mi * mj * (FTable[t] + (ksi - t) * (FTable[t + 1] - FTable[t]));
  return F * rij;
}

void P3MMethod::updateSRForces(int i,
                               int j,
                               int q,
                               int qn,
                               int qnLocal,
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
    shortRangeij = shortRangeForce(rij, particles[i].mass, particles[j].mass);
  } else {
    shortRangeij = shortRangeForceFromTable(rij, particles[i].mass, particles[j].mass);
  }

  particles[i].shortRangeForce += shortRangeij;
  if (qn != q) {
    particles[j].shortRangeFromNeighbor[qnLocal] += -1 * shortRangeij;
  }
}

void P3MMethod::initSRForceTable() {
  for (int i = 0; i < tabulatedValuesCnt; ++i) {
    float rSquared = i * deltaSquared;
    float r = std::sqrtf(rSquared);
    float R;
    if (cloudShape == CloudShape::S1) {
      R = -referenceForceS1(r);
    } else if (cloudShape == CloudShape::S2) {
      R = -referenceForceS2(r);
    } else {
      throw std::invalid_argument("not implemented");
    }

    float G = 1 / (4 * std::numbers::pi_v<float>);
    float eps = softeningLength;
    float totalForce = -G / (r * r + eps * eps);

    FTable.push_back(r == 0 ? 0 : (totalForce - R) / std::sqrtf(r * r + eps * eps));
  }
}

void P3MMethod::updateSRForcesThreadJob(int tid, int threadsCnt, std::vector<Particle>& particles) {
  auto start = std::chrono::steady_clock::now();
  for (int q = tid; q < chainingMesh.getSize(); q += threadsCnt) {
    auto neighbors = chainingMesh.getNeighborsAndSelf(q);
    for (int i = 0; i < neighbors.size(); ++i) {
      int qn = neighbors[i];
      if (qn == -1) {
        continue;
      }

      for (auto node = chainingMesh.getParticlesInCell(q); node != nullptr; node = node->next) {
        for (auto nodeN = chainingMesh.getParticlesInCell(qn); nodeN != nullptr;
             nodeN = nodeN->next) {
          if (enableYSorting &&
              particles[nodeN->particleId].position.y - particles[node->particleId].position.y >
                  cutoffRadius) {
            break;
          }
          updateSRForces(node->particleId, nodeN->particleId, q, qn, i, particles);
        }
      }
    }
  }
  auto end = std::chrono::steady_clock::now();
  auto elapsedMilliSeconds = (end - start).count() * 1e-6f;
  threadTimes[tid] += elapsedMilliSeconds;
}

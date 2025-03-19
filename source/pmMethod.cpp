#include "pmMethod.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <execution>
#include <iostream>
#include <numbers>
#include <ranges>
#include <stdexcept>
#include "greensFunctions.h"
#include "grid.h"
#include "leapfrog.h"
#include "measureTime.h"
#include "simInfo.h"
#include "stateRecorder.h"
#include "unitConversions.h"
#include "vec3.h"

bool isWithingBox(const Vec3& pos, float boxSize) {
  return pos.x >= 0 && pos.x <= boxSize && pos.y >= 0 && pos.y <= boxSize && pos.z >= 0 &&
         pos.z <= boxSize;
}

Vec3 getFieldInCell(int x, int y, int z, FiniteDiffScheme fds, Grid& grid) {
  float fieldX, fieldY, fieldZ;

  if (fds == FiniteDiffScheme::TWO_POINT) {
    fieldX = -0.5f * (grid.getPotential(x + 1, y, z) - grid.getPotential(x - 1, y, z));
    fieldY = -0.5f * (grid.getPotential(x, y + 1, z) - grid.getPotential(x, y - 1, z));
    fieldZ = -0.5f * (grid.getPotential(x, y, z + 1) - grid.getPotential(x, y, z - 1));
  } else if (fds == FiniteDiffScheme::FOUR_POINT) {
    float alpha = 4.0f / 3;
    fieldX = (-1.0f / 12) * (-grid.getPotential(x + 2, y, z) + 8 * grid.getPotential(x + 1, y, z) -
                             8 * grid.getPotential(x - 1, y, z) + grid.getPotential(x - 2, y, z));
    fieldY = (-1.0f / 12) * (-grid.getPotential(x, y + 2, z) + 8 * grid.getPotential(x, y + 1, z) -
                             8 * grid.getPotential(x, y - 1, z) + grid.getPotential(x, y - 2, z));
    fieldZ = (-1.0f / 12) * (-grid.getPotential(x, y, z + 2) + 8 * grid.getPotential(x, y, z + 1) -
                             8 * grid.getPotential(x, y, z - 1) + grid.getPotential(x, y, z - 2));
  } else {
    throw std::invalid_argument("Unknown finite difference type.");
  }

  return Vec3(fieldX, fieldY, fieldZ);
}

void PMMethod::findFourierPotential() {
  int dim = grid.getGridPoints();
  grid.setPotentialFourier(0, 0, 0, std::complex<float>(0, 0));
  auto gridRange = std::ranges::views::iota(0, dim * dim * dim);
  std::for_each(std::execution::par, gridRange.begin(), gridRange.end(), [dim, this](int idx) {
    auto [kx, ky, kz] = grid.indexTripleFromFlat(idx);

    auto densityFourier = grid.getDensityFourier(kx, ky, kz);
    auto potentialFourier = densityFourier * grid.getGreensFunction(kx, ky, kz);
    grid.setPotentialFourier(kx, ky, kz, potentialFourier);
  });
}

void PMMethod::findFieldInCells() {
  auto gridRange = std::ranges::views::iota(0, grid.getLength());
  std::for_each(std::execution::par_unseq, gridRange.begin(), gridRange.end(), [this](int idx) {
    auto [x, y, z] = grid.indexTripleFromFlat(idx);
    auto [fieldX, fieldY, fieldZ] = getFieldInCell(x, y, z, fds, grid);

    Vec3 fieldStrength(fieldX, fieldY, fieldZ);
    grid.assignField(x, y, z, fieldStrength);
  });
}

void PMMethod::updateAccelerations() {
  std::for_each(std::execution::par_unseq, particles.begin(), particles.end(), [this](Particle& p) {
    p.acceleration =
        interpolateField(p.position) +
        accelerationToCodeUnits(externalField(positionToOriginalUnits(p.position, H)), H, DT);
  });
}

void PMMethod::initGreensFunction() {
  int dim = grid.getGridPoints();
  auto gridRange = std::ranges::views::iota(0, dim * dim * dim);
  std::for_each(std::execution::par, gridRange.begin(), gridRange.end(), [dim, this](int idx) {
    auto [kx, ky, kz] = grid.indexTripleFromFlat(idx);

    std::complex<float> G;
    if (gFunc == GreensFunction::DISCRETE_LAPLACIAN) {
      G = GreenDiscreteLaplacian(kx, ky, kz, dim);
    } else if (gFunc == GreensFunction::S1_OPTIMAL) {
      if (is == InterpolationScheme::TSC) {
        G = GreenOptimalTSC(kx, ky, kz, dim, particleDiameter, CloudShape::S1, fds);
      } else {
        throw std::invalid_argument("not implemented");
      }
    } else if (gFunc == GreensFunction::S2_OPTIMAL) {
      if (is == InterpolationScheme::TSC) {
        G = GreenOptimalTSC(kx, ky, kz, dim, particleDiameter, CloudShape::S2, fds);
      } else {
        throw std::invalid_argument("not implemented");
      }
    } else {
      throw std::invalid_argument("not implemented");
    }

    grid.setGreensFunction(kx, ky, kz, G);
  });
}

declareTimeAcc(spreadMass);
declareTimeAcc(forwardFFT);
declareTimeAcc(fourierPotential);
declareTimeAcc(inverseFFT);
declareTimeAcc(fieldInCells);
declareTimeAcc(updateAccelerations);
declareTimeAcc(recordPositions);

void PMMethod::pmMethodStep() {
  measureTime(spreadMass, spreadMass());
  measureTime(forwardFFT, grid.fftDensity());
  measureTime(fourierPotential, findFourierPotential());
  measureTime(inverseFFT, grid.invFftPotential());
  measureTime(fieldInCells, findFieldInCells());
  measureTime(updateAccelerations, updateAccelerations());
}

bool PMMethod::escapedComputationalBox() {
  return std::any_of(
      std::execution::par_unseq, particles.begin(), particles.end(),
      [this](const Particle& p) { return !isWithingBox(p.position, effectiveBoxSize); });
}

Vec3 PMMethod::totalExternalForceOrigUnits() {
  Vec3 f;
  std::for_each(std::execution::seq, particles.begin(), particles.end(),
                [this, &f](const Particle& p) { f += p.mass * externalField(p.position); });

  return f;
}

float TSCAssignmentFunc(float x, int t) {
  if (t == 1) {
    return 0.5f * (0.5f + x) * (0.5f + x);
  }
  if (t == 0) {
    return 0.75f - x * x;
  }
  if (t == -1) {
    return 0.5f * (0.5f - x) * (0.5f - x);
  }
  return 0;
}

void PMMethod::spreadMass() {
  grid.clearDensity();

  switch (is) {
    case InterpolationScheme::NGP: {
      std::for_each(std::execution::par, particles.begin(), particles.end(),
                    [this](const Particle& p) {
                      int x = (int)std::round(p.position.x);
                      int y = (int)std::round(p.position.y);
                      int z = (int)std::round(p.position.z);

                      grid.assignDensity(x, y, z, p.mass);
                    });
      return;
    }

    case InterpolationScheme::CIC: {
      std::for_each(std::execution::par, particles.begin(), particles.end(),
                    [this](const Particle& p) {
                      int x = (int)p.position.x;
                      int y = (int)p.position.y;
                      int z = (int)p.position.z;

                      float d = p.mass;  // divided by vol = H * H * H = 1 (code units)

                      float dx = p.position.x - x;
                      float dy = p.position.y - y;
                      float dz = p.position.z - z;
                      float tx = 1 - dx;
                      float ty = 1 - dy;
                      float tz = 1 - dz;

                      grid.assignDensity(x, y, z, d * tx * ty * tz);
                      grid.assignDensity(x + 1, y, z, d * dx * ty * tz);
                      grid.assignDensity(x, y + 1, z, d * tx * dy * tz);
                      grid.assignDensity(x, y, z + 1, d * tx * ty * dz);

                      grid.assignDensity(x + 1, y + 1, z, d * dx * dy * tz);
                      grid.assignDensity(x + 1, y, z + 1, d * dx * ty * dz);
                      grid.assignDensity(x, y + 1, z + 1, d * tx * dy * dz);

                      grid.assignDensity(x + 1, y + 1, z + 1, d * dx * dy * dz);
                    });

      return;
    }

    case InterpolationScheme::TSC: {
      std::for_each(std::execution::par, particles.begin(), particles.end(),
                    [this](const Particle& p) {
                      int x = (int)p.position.x;
                      int y = (int)p.position.y;
                      int z = (int)p.position.z;

                      float d = p.mass;

                      float dx = p.position.x - x;
                      float dy = p.position.y - y;
                      float dz = p.position.z - z;

                      for (int t1 = -1; t1 <= 1; ++t1) {
                        float T1 = (d / 8) * 2 * TSCAssignmentFunc(dx, t1);
                        for (int t2 = -1; t2 <= 1; ++t2) {
                          float T2 = T1 * 2 * TSCAssignmentFunc(dy, t2);
                          for (int t3 = -1; t3 <= 1; ++t3) {
                            float T3 = T2 * 2 * TSCAssignmentFunc(dz, t3);
                            grid.assignDensity(x + t1, y + t2, z + t3, T3);
                          }
                        }
                      }
                    });
      return;
    }

    default:
      throw std::invalid_argument("Unkown interpolation scheme");
  }
}

Vec3 PMMethod::interpolateField(Vec3 position) {
  float x = position.x;
  float y = position.y;
  float z = position.z;
  switch (is) {
    case InterpolationScheme::NGP: {
      int xi = (int)std::round(x);
      int yi = (int)std::round(y);
      int zi = (int)std::round(z);
      return grid.getField(xi, yi, zi);
    }

    case InterpolationScheme::CIC: {
      int xi = (int)x;
      int yi = (int)y;
      int zi = (int)z;
      float dx = x - xi;
      float dy = y - yi;
      float dz = z - zi;
      float tx = 1 - dx;
      float ty = 1 - dy;
      float tz = 1 - dz;

      return tx * ty * tz * grid.getField(xi, yi, zi) +
             dx * ty * tz * grid.getField(xi + 1, yi, zi) +
             tx * dy * tz * grid.getField(xi, yi + 1, zi) +
             tx * ty * dz * grid.getField(xi, yi, zi + 1) +
             dx * dy * tz * grid.getField(xi + 1, yi + 1, zi) +
             dx * ty * dz * grid.getField(xi + 1, yi, zi + 1) +
             tx * dy * dz * grid.getField(xi, yi + 1, zi + 1) +
             dx * dy * dz * grid.getField(xi + 1, yi + 1, zi + 1);
    }

    case InterpolationScheme::TSC: {
      int xi = (int)x;
      int yi = (int)y;
      int zi = (int)z;
      float dx = x - xi;
      float dy = y - yi;
      float dz = z - zi;

      Vec3 interpolatedField;
      for (int t1 = -1; t1 <= 1; ++t1) {
        float T1 = 2 * TSCAssignmentFunc(dx, t1);
        for (int t2 = -1; t2 <= 1; ++t2) {
          float T2 = T1 * 2 * TSCAssignmentFunc(dy, t2);
          for (int t3 = -1; t3 <= 1; ++t3) {
            float T3 = T2 * 2 * TSCAssignmentFunc(dz, t3);
            interpolatedField += 0.125f * T3 * grid.getField(xi + t1, yi + t2, zi + t3);
          }
        }
      }

      return interpolatedField;
    }

    default:
      throw std::invalid_argument("Unkown interpolation scheme");
  }
}

PMMethod::PMMethod(const std::vector<Vec3>& state,
                   const std::vector<float>& masses,
                   const float effectiveBoxSize,
                   const std::function<Vec3(Vec3)> externalField,
                   const std::function<float(Vec3)> externalPotential,
                   const float H,
                   const float DT,
                   const float G,
                   const InterpolationScheme is,
                   const FiniteDiffScheme fds,
                   const GreensFunction gFunc,
                   const float particleDiameter,
                   Grid& grid)
    : effectiveBoxSize(effectiveBoxSize),
      externalField(externalField),
      externalPotential(externalPotential),
      H(H),
      DT(DT),
      G(G),
      is(is),
      fds(fds),
      gFunc(gFunc),
      particleDiameter(lengthToCodeUnits(particleDiameter, H)),
      grid(grid) {
  this->N = static_cast<int>(masses.size());
  for (int i = 0; i < N; ++i) {
    this->particles.emplace_back(state[i], state[N + i], masses[i]);
  }
}

std::string PMMethod::run(const int simLength,
                          bool collectDiagnostics,
                          bool recordField,
                          const char* positionsPath,
                          const char* energyPath,
                          const char* momentumPath,
                          const char* expectedMomentumPath,
                          const char* fieldPath) {
  StateRecorder stateRecorder(positionsPath, energyPath, momentumPath, expectedMomentumPath,
                              fieldPath);
  SimInfo simInfo;

  if (collectDiagnostics) {
    simInfo.setInitialMomentum(particles);
  }

  stateToCodeUnits(particles, H, DT);
  massToCodeUnits(particles, H, DT, G);
  initGreensFunction();
  pmMethodStep();

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

    measureTime(recordPositions, stateRecorder.recordPositions(particles));
    if (collectDiagnostics) {
      massToOriginalUnits(particles, H, DT, G);
      auto expectedMomentum = simInfo.updateExpectedMomentum(totalExternalForceOrigUnits(), DT);
      stateRecorder.recordExpectedMomentum(expectedMomentum);
      stateRecorder.recordTotalMomentum(SimInfo::totalMomentum(particles));
      auto pe = SimInfo::potentialEnergy(grid, particles, externalPotential, H, DT, G);
      auto ke = SimInfo::kineticEnergy(particles);
      stateRecorder.recordEnergy(pe, ke);
    }
    if (recordField) {
      stateRecorder.recordField(particles, H, DT);
    }
    if (escapedComputationalBox()) {
      std::cout << "Particle moved outside the computational box.\n";
      break;
    }

    stateToCodeUnits(particles, H, DT);
    if (collectDiagnostics) {
      massToCodeUnits(particles, H, DT, G);
      integerStepVelocitiesToCodeUnits(particles, H, DT);
    }

    pmMethodStep();

    updateVelocities(particles);
  }

  printTime(spreadMass);
  printTime(forwardFFT);
  printTime(fourierPotential);
  printTime(inverseFFT);
  printTime(fieldInCells);
  printTime(updateAccelerations);
  printTime(recordPositions);

  return stateRecorder.flush();
}

std::vector<Particle>& PMMethod::getParticles() {
  return particles;
}

float PMMethod::getH() const {
  return H;
}

float PMMethod::getDT() const {
  return DT;
}

float PMMethod::getG() const {
  return G;
}

const Grid& PMMethod::getGrid() const {
  return grid;
}

std::function<float(Vec3)> PMMethod::getExternalPotential() const {
  return externalPotential;
}

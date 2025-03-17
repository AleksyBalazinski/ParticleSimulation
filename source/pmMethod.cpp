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

float sinc(float x) {
  if (x == 0) {
    return 1;
  }
  return std::sinf(x) / x;
}

float TSCFourier(std::array<float, 3> k) {
  float prod = 1;
  for (int i = 0; i < 3; i++) {
    prod *= std::powf(sinc(k[i] / 2), 3);
  }

  return prod;
}

float S1Fourier(float k, float a) {
  float u = k * a / 2;
  return -3 / std::powf(u, 3) * (u * std::cosf(u) - std::sinf(u));
}

float S2Fourier(float k, float a) {
  float u = k * a / 2;
  return 12 / std::powf(u, 4) * (2 - std::cosf(u) - u * std::sinf(u));
}

std::array<std::complex<float>, 3> RFourier(std::array<float, 3> k, float a) {
  std::array<std::complex<float>, 3> R;
  std::complex<float> I(0, 1);
  float kLength = std::sqrtf(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);
  float sSquared = std::powf(S1Fourier(kLength, a), 2);
  for (int i = 0; i < 3; i++) {
    R[i] = -I * float(k[i]) * sSquared / (kLength * kLength);
  }

  return R;
}

std::array<std::complex<float>, 3> DFourier(std::array<float, 3> k) {
  std::array<std::complex<float>, 3> D;
  std::complex<float> I(0, 1);
  for (int i = 0; i < 3; i++) {
    D[i] = I * std::sinf(k[i]);
  }

  return D;
}

std::complex<float> dotProduct(std::array<std::complex<float>, 3> a,
                               std::array<std::complex<float>, 3> b) {
  std::complex<float> dot(0, 0);
  for (int i = 0; i < 3; i++) {
    std::complex<float> biConj(b[i].real(), -1 * b[i].imag());
    dot += a[i] * biConj;
  }

  return dot;
}

void PMMethod::findFourierPotential() {
  int dim = grid.getGridPoints();
  grid.setPotentialFourier(0, 0, 0, std::complex<float>(0, 0));
  auto gridRange = std::ranges::views::iota(0, dim * dim * dim);
  std::for_each(std::execution::par, gridRange.begin(), gridRange.end(), [dim, this](int idx) {
    auto [kx, ky, kz] = grid.indexTripleFromFlat(idx);
    if (kx == 0 && ky == 0 && kz == 0) {
      return;
    }

    std::complex<float> G;
    const float pi = std::numbers::pi_v<float>;
    if (gFunc == GreensFunction::DISCRETE_LAPLACIAN) {
      auto sx = std::sinf(pi * kx / dim);
      auto sy = std::sinf(pi * ky / dim);
      auto sz = std::sinf(pi * kz / dim);
      G = -0.25f / (sx * sx + sy * sy + sz * sz);
    } else if (gFunc == GreensFunction::S1_OPTIMAL) {
      std::array<float, 3> k = {2 * pi * float(kx) / dim, 2 * pi * float(ky) / dim,
                                2 * pi * float(kz) / dim};
      float denomSum = 1;
      for (int i = 0; i < 3; i++) {
        denomSum *=
            (1 - std::powf(std::sinf(k[i] / 2), 2) + 2.0f / 15 * std::powf(std::sinf(k[i] / 2), 4));
      }

      std::array<std::complex<float>, 3> D;
      std::complex<float> I(0, 1);
      if (fds == FiniteDiffScheme::TWO_POINT) {
        D = DFourier(k);
      } else {  // TODO implement four-point
        throw std::invalid_argument("not implemented");
      }

      float DNormSquared = 0;
      for (int i = 0; i < 3; i++) {
        std::complex<float> DiConj(D[i].real(), -1 * D[i].imag());
        DNormSquared += (D[i] * DiConj).real();
      }

      int limit = 2;
      std::array<std::complex<float>, 3> numDotRight = {0, 0, 0};
      for (int n1 = -limit; n1 <= limit; n1++) {
        for (int n2 = -limit; n2 <= limit; n2++) {
          for (int n3 = -limit; n3 <= limit; n3++) {
            std::array<float, 3> kn = {k[0] + 2 * pi * n1, k[1] + 2 * pi * n2, k[2] + 2 * pi * n3};
            float uSquared = std::powf(TSCFourier(kn), 2);

            float a = lengthToCodeUnits(7.5, H);  // TODO
            std::array<std::complex<float>, 3> R = RFourier(kn, a);
            numDotRight[0] += uSquared * R[0];
            numDotRight[1] += uSquared * R[1];
            numDotRight[2] += uSquared * R[2];
          }
        }
      }

      std::complex<float> numerator = dotProduct(D, numDotRight);
      float denominator = DNormSquared * denomSum;

      G = numerator / denominator;
    }

    auto densityFourier = grid.getDensityFourier(kx, ky, kz);
    auto potentialFourier = densityFourier * G;
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

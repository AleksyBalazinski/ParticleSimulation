#include "pmMethod.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <execution>
#include <iostream>
#include <numbers>
#include <ranges>
#include <stdexcept>
#include "grid.h"
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
  std::for_each(std::execution::par_unseq, gridRange.begin(), gridRange.end(),
                [dim, this](int idx) {
                  int kx = idx % dim;
                  int ky = (idx / dim) % dim;
                  int kz = idx / (dim * dim);
                  if (kx == 0 && ky == 0 && kz == 0) {
                    return;
                  }
                  auto sx = std::sinf(std::numbers::pi_v<float> * kx / dim);
                  auto sy = std::sinf(std::numbers::pi_v<float> * ky / dim);
                  auto sz = std::sinf(std::numbers::pi_v<float> * kz / dim);
                  auto G = -0.25f / (sx * sx + sy * sy + sz * sz);

                  auto densityFourier = grid.getDensityFourier(kx, ky, kz);
                  auto potentialFourier =
                      std::complex<float>(G * densityFourier.real(), G * densityFourier.imag());
                  grid.setPotentialFourier(kx, ky, kz, potentialFourier);
                });
}

void PMMethod::findFieldInCells() {
  int dim = grid.getGridPoints();
  auto gridRange = std::ranges::views::iota(0, dim * dim * dim);
  std::for_each(std::execution::par_unseq, gridRange.begin(), gridRange.end(),
                [dim, this](int idx) {
                  int x = idx % dim;
                  int y = (idx / dim) % dim;
                  int z = idx / (dim * dim);
                  auto [fieldX, fieldY, fieldZ] = getFieldInCell(x, y, z, fds, grid);

                  Vec3 fieldStrength(fieldX, fieldY, fieldZ);
                  grid.assignField(x, y, z, fieldStrength);
                });
}

void PMMethod::updateAccelerations() {
  auto nRange = std::ranges::views::iota(0, N);
  std::for_each(std::execution::par_unseq, nRange.begin(), nRange.end(), [this](int i) {
    accelerations[i] =
        interpolateField(state[i].x, state[i].y, state[i].z) + externalField(state[i]);
  });
}

void PMMethod::pmMethodStep() {
  reassignDensity();
  grid.fftDensity();
  findFourierPotential();
  grid.invFftPotential();
  findFieldInCells();
  updateAccelerations();
}

bool PMMethod::escapedComputationalBox() {
  return std::any_of(std::execution::par_unseq, state.begin(), state.begin() + N,
                     [this](const Vec3 pos) { return !isWithingBox(pos, effectiveBoxSize); });
}

void PMMethod::reassignDensity() {
  auto nIdxRange = std::ranges::views::iota(0, N);
  grid.clearDensity();

  switch (is) {
    case InterpolationScheme::NGP:
      std::for_each(std::execution::par_unseq, nIdxRange.begin(), nIdxRange.end(), [this](int i) {
        int x = (int)std::round(state[i].x);
        int y = (int)std::round(state[i].y);
        int z = (int)std::round(state[i].z);

        float vol = H * H * H;
        grid.assignDensity(x, y, z, densityToCodeUnits(masses[i] / vol, DT, G));
      });
      return;

    case InterpolationScheme::CIC:
      std::for_each(std::execution::unseq, nIdxRange.begin(), nIdxRange.end(), [this](int i) {
        int x = (int)state[i].x;
        int y = (int)state[i].y;
        int z = (int)state[i].z;

        float vol = H * H * H;
        float d = densityToCodeUnits(masses[i] / vol, DT, G);

        float dx = state[i].x - x;
        float dy = state[i].y - y;
        float dz = state[i].z - z;
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
}

Vec3 PMMethod::interpolateField(float x, float y, float z) {
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

    default:
      throw std::invalid_argument("Unkown interpolation scheme");
  }
}

void setIntegerVelocities(std::vector<Vec3>& intVs,
                          const std::vector<Vec3>& state,
                          const std::vector<float>& masses,
                          float G,
                          float h,
                          const std::vector<Vec3>& accelerations) {
  int n = (int)intVs.size();
  for (int i = 0; i < n; i++) {
    intVs[i] = state[n + i] + 0.5f * h * accelerations[i];
  }
}

PMMethod::PMMethod(std::vector<Vec3>& state,
                   std::vector<float>& masses,
                   float effectiveBoxSize,
                   std::function<Vec3(Vec3)> externalField,
                   float H,
                   float DT,
                   float G,
                   InterpolationScheme is,
                   FiniteDiffScheme fds,
                   Grid& grid)
    : state(state),
      masses(masses),
      effectiveBoxSize(effectiveBoxSize),
      H(H),
      DT(DT),
      G(G),
      is(is),
      fds(fds),
      grid(grid) {
  this->externalField = [DT, H, externalField](Vec3 pos) -> Vec3 {
    return accelerationToCodeUnits(externalField(positionToOriginalUnits(pos, H)), H, DT);
  };
  this->N = static_cast<int>(masses.size());
  this->accelerations.resize(N);
  this->intStepVelocities.resize(N);
}

void PMMethod::setHalfVelocities() {
  auto nRange = std::ranges::views::iota(0, N);
  std::for_each(std::execution::par_unseq, nRange.begin(), nRange.end(),
                [this](int i) { state[N + i] += 0.5 * accelerations[i]; });
}

void PMMethod::updateVelocities() {
  auto nRange = std::ranges::views::iota(0, N);
  std::for_each(std::execution::par_unseq, nRange.begin(), nRange.end(),
                [this](int i) { state[N + i] += accelerations[i]; });
}

void PMMethod::updatePositions() {
  auto nRange = std::ranges::views::iota(0, N);
  std::for_each(std::execution::par_unseq, nRange.begin(), nRange.end(),
                [this](int i) { state[i] += state[N + i]; });
}

std::string PMMethod::run(const int simLength,
                          bool collectDiagnostics,
                          const char* positionsPath,
                          const char* energyPath,
                          const char* momentumPath) {
  StateRecorder stateRecorder(positionsPath, energyPath, momentumPath);

  stateToCodeUnits(state, H, DT);
  pmMethodStep();

  setHalfVelocities();

  for (int t = 0; t <= simLength; ++t) {
    std::cout << "progress: " << float(t) / simLength << '\r';
    std::cout.flush();

    updatePositions();

    // if (collectDiagnostics) {
    //   setIntegerVelocities(intStepVelocities, state, masses, G, 1, accelerations);
    // }
    stateToOriginalUnits(state, H, DT);
    // if (collectDiagnostics) {
    //   velocitiesToOriginalUnits(intStepVelocities, H, DT);
    // }

    stateRecorder.recordPositions(state.begin(), state.begin() + N);
    // if (collectDiagnostics) {
    //   stateRecorder.recordEnergy(
    //       potentialEnergy(state.begin(), state.begin() + N, masses, G),
    //       kineticEnergy(intStepVelocities.begin(), intStepVelocities.end(), masses, G));

    //   stateRecorder.recordTotalMomentum(
    //       totalMomentum(intStepVelocities.begin(), intStepVelocities.end(), masses));
    // }
    if (escapedComputationalBox()) {
      std::cout << "Particle moved outside the computational box.\n";
      break;
    }

    stateToCodeUnits(state, H, DT);
    // if (collectDiagnostics) {
    //   velocitiesToCodeUnits(intStepVelocities, H, DT);
    // }

    pmMethodStep();

    updateVelocities();
  }
  return stateRecorder.flush();
}
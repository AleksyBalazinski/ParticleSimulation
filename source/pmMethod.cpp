#include "pmMethod.h"
#include <kiss_fftnd.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <execution>
#include <iostream>
#include <ranges>
#include <stdexcept>
#include "grid.h"
#include "simInfo.h"
#include "stateRecorder.h"
#include "unit_conversions.h"
#include "vec3.h"

// #define DEBUG

#ifdef DEBUG
long long totalTimeMs = 0;
long long fftTimeMs = 0;
#endif

bool isWithingBox(const Vec3& pos, double boxSize) {
  return pos.x >= 0 && pos.x <= boxSize && pos.y >= 0 && pos.y <= boxSize && pos.z >= 0 &&
         pos.z <= boxSize;
}

void PMMethod::updateAccelerations(std::vector<Vec3>& accelerations,
                                   const std::vector<Vec3>& state,
                                   const std::vector<double>& masses,
                                   double H,
                                   double DT,
                                   double G) {
#ifdef DEBUG
  auto beginAll = std::chrono::steady_clock::now();
#endif
  int n = (int)masses.size();
  if (int boxSize = grid.getGridPoints();
      std::any_of(state.begin(), state.begin() + n,
                  [boxSize](const Vec3 pos) { return !isWithingBox(pos, boxSize); })) {
    throw std::runtime_error("A particle moved outside the computational box.");
  }

  reassignDensity(state, masses, H, DT, G);

#ifdef DEBUG
  auto beginFFT1 = std::chrono::steady_clock::now();
#endif
  grid.fftDensity();
#ifdef DEBUG
  auto endFFT1 = std::chrono::steady_clock::now();
#endif

  // find potential in Fourier space
  int dim = grid.getGridPoints();
  grid.setPotentialFourier(0, 0, 0, kiss_fft_cpx(0, 0));
  auto gridIdxRange = std::ranges::views::iota(0, dim * dim * dim);
  std::for_each(std::execution::par, gridIdxRange.begin(), gridIdxRange.end(),
                [dim, this](int idx) {
                  int kx = idx / (dim * dim);
                  int ky = (idx / dim) % dim;
                  int kz = idx % dim;
                  if (kx == 0 && ky == 0 && kz == 0) {
                    return;
                  }
                  auto sx = std::sin(std::numbers::pi * kx / dim);
                  auto sy = std::sin(std::numbers::pi * ky / dim);
                  auto sz = std::sin(std::numbers::pi * kz / dim);
                  double G = -0.25 / (sx * sx + sy * sy + sz * sz);

                  auto densityFourier = grid.getDensityFourier(kx, ky, kz);
                  auto potentialFourier = kiss_fft_cpx(G * densityFourier.r, G * densityFourier.i);
                  grid.setPotentialFourier(kx, ky, kz, potentialFourier);
                });

  // find real potential by applying IFFT
#ifdef DEBUG
  auto beginFFT2 = std::chrono::steady_clock::now();
#endif
  grid.invFftPotential();

#ifdef DEBUG
  auto endFFT2 = std::chrono::steady_clock::now();
#endif

  // find field in a meshpoint
  std::for_each(
      std::execution::par, gridIdxRange.begin(), gridIdxRange.end(), [dim, this](int idx) {
        int x = idx / (dim * dim);
        int y = (idx / dim) % dim;
        int z = idx % dim;
        auto fieldX = -0.5 * (grid.getPotential(x + 1, y, z) - grid.getPotential(x - 1, y, z));
        auto fieldY = -0.5 * (grid.getPotential(x, y + 1, z) - grid.getPotential(x, y - 1, z));
        auto fieldZ = -0.5 * (grid.getPotential(x, y, z + 1) - grid.getPotential(x, y, z - 1));

        Vec3 fieldStrength(fieldX, fieldY, fieldZ);
        grid.assignField(x, y, z, fieldStrength);
      });

  // acceleration calculation
  for (int i = 0; i < n; ++i) {
    accelerations[i] = getFieldAtMeshpoint(state[i].x, state[i].y, state[i].z);
  }
#ifdef DEBUG
  auto endAll = std::chrono::steady_clock::now();

  totalTimeMs += std::chrono::duration_cast<std::chrono::milliseconds>(endAll - beginAll).count();
  fftTimeMs += std::chrono::duration_cast<std::chrono::milliseconds>(endFFT1 - beginFFT1).count() +
               std::chrono::duration_cast<std::chrono::milliseconds>(endFFT2 - beginFFT2).count();
#endif
}

void PMMethod::reassignDensity(const std::vector<Vec3>& state,
                               const std::vector<double>& masses,
                               double H,
                               double DT,
                               double G) {
  int n = (int)masses.size();
  grid.clearDensity();
  for (int i = 0; i < n; i++) {
    int x = (int)std::round(state[i].x);
    int y = (int)std::round(state[i].y);
    int z = (int)std::round(state[i].z);
    double vol = H * H * H;
    grid.assignDensity(x, y, z, densityToCodeUnits(masses[i] / vol, DT, G));
  }
}

Vec3 PMMethod::getFieldAtMeshpoint(double x, double y, double z) {
  int xi = (int)std::round(x);
  int yi = (int)std::round(y);
  int zi = (int)std::round(z);

  return grid.getField(xi, yi, zi);
}

void setIntegerVelocities(std::vector<Vec3>& intVs,
                          const std::vector<Vec3>& state,
                          const std::vector<double>& masses,
                          double G,
                          double h,
                          const std::vector<Vec3>& accelerations) {
  int n = (int)intVs.size();
  for (int i = 0; i < n; i++) {
    intVs[i] = state[n + i] + 0.5 * h * accelerations[i];
  }
}

PMMethod::PMMethod(int gridPoints) : grid(gridPoints) {}

std::string PMMethod::run(std::vector<Vec3>& state,
                          std::vector<double>& masses,
                          const double simLengthSeconds,
                          const double stepSize,
                          const double cellSize,
                          const double G,
                          const int frameRate,
                          const char* outPath,
                          const char* energyPath,
                          const char* momentumPath) {
  StateRecorder stateRecorder(outPath, energyPath, momentumPath);
  const int n = (int)masses.size();
  double curFrameAcc = 0;
  const double frameLength = 1.0 / frameRate;
  std::vector<Vec3> accelerations(n);
  std::vector<Vec3> velocities(n);  // velocities at integer step (needed only for display)

  stateToCodeUnits(state, cellSize, stepSize);
  updateAccelerations(accelerations, state, masses, cellSize, stepSize, G);

  // set v_(1/2)
  // from this point on state[n + i] holds velocities at half-step
  for (int i = 0; i < n; i++) {
    state[n + i] += 0.5 * accelerations[i];
  }

  for (double t = 0; t <= simLengthSeconds; t += stepSize) {
    std::cout << "progress: " << t / simLengthSeconds << '\r';
    std::cout.flush();
    for (int i = 0; i < n; i++) {
      state[i] += state[n + i];
    }
    if (curFrameAcc <= 0) {
      setIntegerVelocities(velocities, state, masses, G, 1, accelerations);
      stateToOriginalUnits(state, cellSize, stepSize);
      velocitiesToOriginalUnits(velocities, cellSize, stepSize);
      stateRecorder.recordPositions(state.begin(), state.begin() + n);
      stateRecorder.recordEnergy(potentialEnergy(state.begin(), state.begin() + n, masses, G),
                                 kineticEnergy(velocities.begin(), velocities.end(), masses, G));

      stateRecorder.recordTotalMomentum(
          totalMomentum(velocities.begin(), velocities.end(), masses));
      curFrameAcc = frameLength;
      stateToCodeUnits(state, cellSize, stepSize);
      velocitiesToCodeUnits(velocities, cellSize, stepSize);
    }
    updateAccelerations(accelerations, state, masses, cellSize, stepSize, G);
    // now that we have accelerations of all particles, we can predict motion
    for (int i = 0; i < n; i++) {
      state[n + i] += accelerations[i];
    }
    curFrameAcc -= stepSize;
  }

#ifdef DEBUG
  std::cout << "total time: " << totalTimeMs << "[ms]\n";
  std::cout << "FFT time: " << fftTimeMs << "[ms]\n";
  std::cout << "FFT contrib.: " << (float)fftTimeMs / totalTimeMs * 100 << "%\n";
#endif
  return stateRecorder.flush();
}
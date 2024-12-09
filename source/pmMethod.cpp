#include "pmMethod.h"
#include <kiss_fftnd.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <execution>
#include <iostream>
#include <numbers>
#include <ranges>
#include "grid.h"
#include "simInfo.h"
#include "stateRecorder.h"
#include "vec3.h"

#define DEBUG

#ifdef DEBUG
long long totalTimeMs = 0;
long long fftTimeMs = 0;
#endif

void PMMethod::updateAccelerations(double G) {
#ifdef DEBUG
  auto beginAll = std::chrono::steady_clock::now();
#endif
  int n = (int)masses.size();
  reassignDensity(masses, G);

#ifdef DEBUG
  auto beginFFT1 = std::chrono::steady_clock::now();
#endif
  auto densityFourier = grid.fftDensity();
#ifdef DEBUG
  auto endFFT1 = std::chrono::steady_clock::now();
#endif

  // find potential in Fourier space
  int dim = grid.getGridPoints();
  grid.setPotentialFourier(0, 0, 0, kiss_fft_cpx(0, 0));
  auto gridIdxRange = std::ranges::views::iota(0, dim * dim * dim);
  std::for_each(std::execution::par, gridIdxRange.begin(), gridIdxRange.end(), [&](int idx) {
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
  auto potential = grid.invFftPotential();

#ifdef DEBUG
  auto endFFT2 = std::chrono::steady_clock::now();
#endif

  // find field in a meshpoint
  std::for_each(std::execution::par, gridIdxRange.begin(), gridIdxRange.end(), [&](int idx) {
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

PMMethod::PMMethod(std::vector<Vec3>& state,
                   std::vector<double>& masses,
                   int gridPoints,
                   double H,
                   double DT)
    : state(state),
      masses(masses),
      accelerations(masses.size()),
      velocities(masses.size()),
      gridPoints(gridPoints),
      N(masses.size()),
      H(H),
      DT(DT),
      grid(gridPoints) {}

Vec3 PMMethod::positionInCodeUntits(const Vec3& pos) {
  return pos / H;
}

Vec3 PMMethod::velocityInCodeUnits(const Vec3& v) {
  return DT * v / H;
}

double PMMethod::densityToCodeUnits(double density, double G) {
  return DT * DT * 4 * std::numbers::pi * G * density;
}

void PMMethod::stateToCodeUnits() {
  std::transform(state.begin(), state.begin() + N, state.begin(),
                 [this](const Vec3& pos) { return positionInCodeUntits(pos); });
  std::transform(state.begin() + N, state.end(), state.begin() + N,
                 [this](const Vec3& v) { return velocityInCodeUnits(v); });
}

void PMMethod::velocitiesToCodeUnits() {
  std::transform(velocities.begin(), velocities.end(), state.begin() + N,
                 [this](const Vec3& v) { return velocityInCodeUnits(v); });
}

Vec3 PMMethod::positionInOriginalUnits(const Vec3& pos) {
  return H * pos;
}

Vec3 PMMethod::velocityInOriginalUnits(const Vec3& v) {
  return H * v / DT;
}

void PMMethod::stateToOriginalUnits() {
  std::transform(state.begin(), state.begin() + N, state.begin(),
                 [this](const Vec3& pos) { return positionInOriginalUnits(pos); });
  std::transform(state.begin() + N, state.end(), state.begin() + N,
                 [this](const Vec3& v) { return velocityInOriginalUnits(v); });
}

void PMMethod::velocitiesToOriginalUnits() {
  std::transform(velocities.begin(), velocities.end(), velocities.begin(),
                 [this](const Vec3& v) { return velocityInOriginalUnits(v); });
}

void PMMethod::reassignDensity(const std::vector<double>& masses, double G) {
  int n = (int)masses.size();
  grid.clearDensity();
  for (int i = 0; i < n; i++) {
    int x = (int)std::round(state[i].x);
    int y = (int)std::round(state[i].y);
    int z = (int)std::round(state[i].z);
    grid.assignDensity(
        x, y, z,
        densityToCodeUnits(masses[i], G));  // account for units change TODO: divide by volume??
  }
}

Vec3 PMMethod::getFieldAtMeshpoint(double x, double y, double z) {
  int xi = (int)std::round(x);
  int yi = (int)std::round(y);
  int zi = (int)std::round(z);

  return grid.getField(xi, yi, zi);
}

void setIntegerVelocities(const std::vector<Vec3>& state,
                          const std::vector<double>& masses,
                          double G,
                          double h,
                          const std::vector<Vec3>& accelerations,
                          std::vector<Vec3>& intVs) {
  int n = (int)intVs.size();
  for (size_t i = 0; i < n; i++) {
    intVs[i] = state[n + i] + 0.5 * h * accelerations[i];
  }
}

std::string PMMethod::run(const double simLengthSeconds,
                          const double G,
                          const int frameRate,
                          const char* outPath,
                          const char* energyPath,
                          const char* momentumPath) {
  StateRecorder stateRecorder(outPath, energyPath, momentumPath);
  const int n = (int)masses.size();
  double curFrameAcc = 0;
  const double frameLength = 1.0 / frameRate;

  stateToCodeUnits();
  updateAccelerations(G);

  // set v_(1/2)
  // from this point on state[n + i] holds velocities at half-step
  for (int i = 0; i < n; i++) {
    state[n + i] += 0.5 * accelerations[i];
  }

  for (double t = 0; t <= simLengthSeconds; t += DT) {
    std::cout << "progress: " << t / simLengthSeconds << '\r';
    std::cout.flush();
    for (int i = 0; i < n; i++) {
      state[i] += state[n + i];
    }
    if (curFrameAcc <= 0) {
      setIntegerVelocities(state, masses, G, 1, accelerations, velocities);
      stateToOriginalUnits();
      velocitiesToOriginalUnits();
      stateRecorder.recordState(state.begin(), state.begin() + n);
      stateRecorder.recordTotalEnergy(totalEnergy(state.begin(), state.begin() + n,
                                                  velocities.begin(), velocities.end(), masses, G));

      stateRecorder.recordTotalMomentum(
          totalMomentum(velocities.begin(), velocities.end(), masses));
      curFrameAcc = frameLength;
      stateToCodeUnits();
      velocitiesToCodeUnits();
    }
    updateAccelerations(G);
    // now that we have accelerations of all particles, we can predict motion
    for (int i = 0; i < n; i++) {
      state[n + i] += accelerations[i];
    }
    // TODO: what if a particle ends up outside the grid?
    curFrameAcc -= DT;
  }

#ifdef DEBUG
  std::cout << "total time: " << totalTimeMs << "[ms]\n";
  std::cout << "FFT time: " << fftTimeMs << "[ms]\n";
  std::cout << "FFT contrib.: " << (float)fftTimeMs / totalTimeMs * 100 << "%\n";
#endif
  return stateRecorder.flush();
}
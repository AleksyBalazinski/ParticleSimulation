#include "pmMethod.h"
#include <kiss_fftnd.h>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numbers>
#include "grid.h"
#include "simInfo.h"
#include "stateRecorder.h"
#include "vec3.h"

// #define DEBUG

void PMMethod::updateAccelerations(double G) {
#ifdef DEBUG
  auto beginAll = std::chrono::steady_clock::now();
#endif
  int n = (int)masses.size();
  int nfft = grid.getLength();
  reassignDensity(masses, G);

  for (int i = 0; i < nfft; ++i) {
    density[i].r = grid.getDensity()[i];
  }

#ifdef DEBUG
  auto beginFFT1 = std::chrono::steady_clock::now();
#endif
  kiss_fftnd(cfg, density, densityFourier);
#ifdef DEBUG
  auto endFFT1 = std::chrono::steady_clock::now();
#endif

  // find potential in Fourier space
  int dim = grid.getGridPoints();
  for (int kx = 0; kx < dim; ++kx) {
    for (int ky = 0; ky < dim; ++ky) {
      for (int kz = 0; kz < dim; ++kz) {
        if (kx == 0 && ky == 0 && kz == 0) {
          potentialFourier[0].r = 0;
          potentialFourier[0].i = 0;
          continue;
        }
        auto sx = std::sin(std::numbers::pi * kx / dim);
        auto sy = std::sin(std::numbers::pi * ky / dim);
        auto sz = std::sin(std::numbers::pi * kz / dim);
        double G = -0.25 / (sx * sx + sy * sy + sz * sz);

        int idx = kx * dim * dim + ky * dim + kz;
        potentialFourier[idx].r = G * densityFourier[idx].r;
        potentialFourier[idx].i = G * densityFourier[idx].i;
      }
    }
  }

  // find real potential by applying IFFT
#ifdef DEBUG
  auto beginFFT2 = std::chrono::steady_clock::now();
#endif
  kiss_fftnd(cfgInv, potentialFourier, potential);
  for (int i = 0; i < nfft; i++) {
    potential[i].r /= nfft;
    potential[i].i /= nfft;
  }

#ifdef DEBUG
  auto endFFT2 = std::chrono::steady_clock::now();
#endif

  // find field in a meshpoint
  for (int x = 0; x < dim; ++x) {
    for (int y = 0; y < dim; ++y) {
      for (int z = 0; z < dim; ++z) {
        auto fieldX = -0.5 * (potential[getIndx(x + 1, y, z, dim)].r -
                              potential[getIndx(x - 1, y, z, dim)].r);
        auto fieldY = -0.5 * (potential[getIndx(x, y + 1, z, dim)].r -
                              potential[getIndx(x, y - 1, z, dim)].r);
        auto fieldZ = -0.5 * (potential[getIndx(x, y, z + 1, dim)].r -
                              potential[getIndx(x, y, z - 1, dim)].r);

        Vec3 fieldStrength(fieldX, fieldY, fieldZ);
        grid.assignField(x, y, z, fieldStrength);
      }
    }
  }

  // acceleration calculation
  for (int i = 0; i < n; ++i) {
    accelerations[i] = getFieldAtMeshpoint(state[i].x, state[i].y, state[i].z);
  }
#ifdef DEBUG
  auto endAll = std::chrono::steady_clock::now();

  std::cout << "Total time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(endAll - beginAll).count()
            << "[ms]\n";
  std::cout
      << "FFT time: "
      << std::chrono::duration_cast<std::chrono::milliseconds>(endFFT1 - beginFFT1).count() +
             std::chrono::duration_cast<std::chrono::milliseconds>(endFFT2 - beginFFT2).count()
      << "[ms]\n";
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
      gridPoints(gridPoints),
      H(H),
      DT(DT),
      grid(gridPoints) {
  int nfft = grid.getLength();

  density = new kiss_fft_cpx[nfft];
  densityFourier = new kiss_fft_cpx[nfft];

  potentialFourier = new kiss_fft_cpx[nfft];
  potential = new kiss_fft_cpx[nfft];

  int dim = grid.getGridPoints();
  int dims[] = {dim, dim, dim};
  cfg = kiss_fftnd_alloc(dims, 3, false, nullptr, nullptr);
  cfgInv = kiss_fftnd_alloc(dims, 3, true, nullptr, nullptr);
}

Vec3 PMMethod::positionInCodeUntits(const Vec3& pos) {
  return pos / H;
}

Vec3 PMMethod::velocityInCodeUnits(const Vec3& v) {
  return DT * v / H;
}

double PMMethod::densityToCodeUnits(double density, double G) {
  return DT * DT * 4 * std::numbers::pi * G * density;
}

void PMMethod::stateToCodeUnits(std::vector<Vec3>& state, int n) {
  for (int i = 0; i < n; i++) {
    state[i] = positionInCodeUntits(state[i]);
    state[n + i] = velocityInCodeUnits(state[n + i]);
  }
}

Vec3 PMMethod::positionInOriginalUnits(const Vec3& pos) {
  return H * pos;
}

Vec3 PMMethod::velocityInOriginalUnits(const Vec3& v) {
  return H * v / DT;
}

void PMMethod::stateToOriginalUnits(std::vector<Vec3>& state, int n) {
  for (int i = 0; i < n; ++i) {
    state[i] = positionInOriginalUnits(state[i]);
    state[n + i] = velocityInOriginalUnits(state[n + i]);
  }
}

void PMMethod::reassignDensity(const std::vector<double>& masses, double G) {
  int n = (int)masses.size();
  grid.clearDensity();
  for (int i = 0; i < n; i++) {
    int x = (int)std::round(state[i].x);
    int y = (int)std::round(state[i].y);
    int z = (int)std::round(state[i].z);
    grid.assignDensity(x, y, z, densityToCodeUnits(masses[i], G));  // account for units change
  }
}

Vec3 PMMethod::getFieldAtMeshpoint(double x, double y, double z) {
  int xi = (int)std::round(x);
  int yi = (int)std::round(y);
  int zi = (int)std::round(z);

  return grid.getField(xi, yi, zi);
}

PMMethod::~PMMethod() {
  kiss_fft_free(cfg);
  delete[] density;
  delete[] densityFourier;

  kiss_fft_free(cfgInv);
  delete[] potential;
  delete[] potentialFourier;
}

void printState(const std::vector<Vec3>& state) {
  for (const auto& i : state) {
    std::cout << i.toString() << '\n';
  }
  std::cout << '\n';
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

  // set v_(1/2)
  // from this point on state[n + i] holds velocities at half-step
  for (int i = 0; i < n; i++) {
    state[n + i] += 0.5 * accelerations[i];
  }
  std::vector<Vec3> velocities(n);  // velocities at integer step (needed only for display)

  stateToCodeUnits(state, n);
  for (double t = 0; t <= simLengthSeconds; t += DT) {
    std::cout << "progress: " << t / simLengthSeconds << '\r';
    std::cout.flush();
    for (int i = 0; i < n; i++) {
      state[i] += state[n + i];
    }
    if (curFrameAcc <= 0) {
      setIntegerVelocities(state, masses, G, 1, accelerations, velocities);
      stateToOriginalUnits(state, n);
      // TODO: velocities should also be converted but here DT = H, so it wouldn't change anything
      stateRecorder.recordState(state.begin(), state.begin() + n);
      stateRecorder.recordTotalEnergy(totalEnergy(state.begin(), state.begin() + n,
                                                  velocities.begin(), velocities.end(), masses, G));

      stateRecorder.recordTotalMomentum(
          totalMomentum(velocities.begin(), velocities.end(), masses));
      curFrameAcc = frameLength;
      stateToCodeUnits(state, n);
      // TODD: velocties should be converted back to code units
    }
    updateAccelerations(G);
    // now that we have accelerations of all particles, we can predict motion
    for (int i = 0; i < n; i++) {
      state[n + i] += accelerations[i];
    }
    // TODO: what if a particle ends up outside the grid?
    curFrameAcc -= DT;
  }

  return stateRecorder.flush();
}
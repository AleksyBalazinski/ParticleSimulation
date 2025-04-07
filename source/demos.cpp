#include "FFTWAdapter.h"
#ifdef CUDA
#include "PMMethodGPU.h"
#endif
#include "diskSamplerLinear.h"
#include "externalFields.h"
#include "grid.h"
#include "kissFFTAdapter.h"
#include "p3mMethod.h"
#include "pmMethod.h"
#include "ppMethod.h"
#include "utils.h"

std::vector<Vec3> linearSpacing(int n, float dx, Vec3 center) {
  std::vector<Vec3> state(2 * n);

  for (int i = 0; i < n; ++i) {
    state[i] = center + Vec3(i * dx, 0, 0);
    state[n + i] = Vec3::zero();
  }

  return state;
}

void probeField() {
  /*
  units:
  [M] = G(solar mass)
  [L] = kpc
  [T] = Myr
  */
  const int n = 400;
  const float G = 4.5e-3f;

  std::vector<Vec3> state = linearSpacing(n, 0.075f /*dx*/, Vec3(30, 30, 15) /*center*/);

  std::vector<float> masses(n, 0);
  masses[0] = 1.0f;

  auto externalField = [](Vec3 pos) -> Vec3 { return Vec3::zero(); };
  auto externalPotential = [](Vec3 pos) -> float { return 0; };

  auto gridPoints = std::make_tuple(64, 64, 32);
  std::array<int, 3> dims = {std::get<2>(gridPoints), std::get<1>(gridPoints),
                             std::get<0>(gridPoints)};
  FFTWAdapter fftAdapter(dims);
  Grid grid(gridPoints, fftAdapter);
  auto effectiveBoxSize = std::make_tuple(60.0f, 60.0f, 30.0f);
  float H = std::get<0>(effectiveBoxSize) / (std::get<0>(gridPoints) / 2);
  float DT = 1;
  float particleDiameter = 4 * H;

  PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
              InterpolationScheme::TSC, FiniteDiffScheme::TWO_POINT, GreensFunction::S1_OPTIMAL,
              particleDiameter, grid);

  pm.run(2, true /*diagnostics*/, true /*field*/);
}

void galaxySimulationPM(int n, int simLength) {
  /*
    units:
    [M] = G(solar mass)
    [L] = kpc
    [T] = Myr
    */
  Vec3 galaxyCenter = Vec3(30, 30, 15);
  float rb = 3.0f;
  float mb = 60.0f;
  float rd = 15.0f;
  float md = 15.0f;
  float thickness = 0.3f;
  float G = 4.5e-3f;

  DiskSamplerLinear diskSampler;
  std::vector<Vec3> state = diskSampler.sample(galaxyCenter, rb, mb, rd, md, thickness, G, n);

  std::vector<float> masses(n, md / n);
  auto externalField = [galaxyCenter, rb, mb, G](Vec3 pos) -> Vec3 {
    return sphRadDecrField(pos, galaxyCenter, rb, mb, G);
  };
  auto externalPotential = [galaxyCenter, rb, mb, G](Vec3 pos) -> float {
    return sphRadDecrFieldPotential(pos, galaxyCenter, rb, mb, G);
  };

  auto gridPoints = std::make_tuple(128, 128, 64);
  std::array<int, 3> dims = {std::get<2>(gridPoints), std::get<1>(gridPoints),
                             std::get<0>(gridPoints)};
  auto effectiveBoxSize = std::make_tuple(60.0f, 60.0f, 30.0f);
  float H = std::get<0>(effectiveBoxSize) / (std::get<0>(gridPoints) / 2);
  float DT = 1;

#ifdef CUDA
  PMMethodGPU pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
                 InterpolationScheme::TSC, FiniteDiffScheme::TWO_POINT,
                 GreensFunction::DISCRETE_LAPLACIAN, 0, gridPoints);
#else
  FFTWAdapter fftAdapter(dims);
  Grid grid(gridPoints, fftAdapter);
  PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
              InterpolationScheme::TSC, FiniteDiffScheme::TWO_POINT,
              GreensFunction::DISCRETE_LAPLACIAN, 0, grid);
#endif
  pm.run(simLength, true /*diagnostics*/);
}

void galaxySimulationP3M(int n, int simLength) {
  /*
    units:
    [M] = G(solar mass)
    [L] = kpc
    [T] = Myr
    */
  Vec3 galaxyCenter = Vec3(30, 30, 15);
  float rb = 3.0f;
  float mb = 60.0f;
  float rd = 15.0f;
  float md = 15.0f;
  float thickness = 0.3f;
  float G = 4.5e-3f;

  DiskSamplerLinear diskSampler;
  std::vector<Vec3> state = diskSampler.sample(galaxyCenter, rb, mb, rd, md, thickness, G, n);

  std::vector<float> masses(n, md / n);
  auto externalField = [galaxyCenter, rb, mb, G](Vec3 pos) -> Vec3 {
    return sphRadDecrField(pos, galaxyCenter, rb, mb, G);
  };
  auto externalPotential = [galaxyCenter, rb, mb, G](Vec3 pos) -> float {
    return sphRadDecrFieldPotential(pos, galaxyCenter, rb, mb, G);
  };

  auto gridPoints = std::make_tuple(128, 128, 64);
  std::array<int, 3> dims = {std::get<2>(gridPoints), std::get<1>(gridPoints),
                             std::get<0>(gridPoints)};
  FFTWAdapter fftAdapter(dims);
  Grid grid(gridPoints, fftAdapter);
  auto effectiveBoxSize = std::make_tuple(60.0f, 60.0f, 30.0f);
  float H = std::get<0>(effectiveBoxSize) / (std::get<0>(gridPoints) / 2);
  float DT = 1;
  float a = 3 * H;

#ifdef CUDA
  PMMethodGPU pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
                 InterpolationScheme::TSC, FiniteDiffScheme::TWO_POINT,
                 GreensFunction::DISCRETE_LAPLACIAN, 0, gridPoints);
#else
  PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
              InterpolationScheme::TSC, FiniteDiffScheme::TWO_POINT, GreensFunction::S1_OPTIMAL, a,
              grid);
#endif

  float re = 0.7f * a;
  float softeningLength = 1.5f;
  P3MMethod p3m(pm, effectiveBoxSize, re, a, H, softeningLength, CloudShape::S1);

  p3m.run(simLength, true /*diagnostics*/);
}

void smallSimPP() {
  const int n = 3;
  std::vector<float> masses = {20, 5, 1e2};
  std::vector<Vec3> state = {
      Vec3(30, 30, 15), Vec3(45, 32, 15),  Vec3(30, 10, 15),  // positions
      Vec3(0.1f, 0, 0), Vec3(-0.3f, 0, 0), Vec3::zero()       // velocities
  };

  int simLength = 100;
  float DT = 1;
  float G = 4.5e-3f;

  ppMethodLeapfrog(state, masses, simLength, DT, G, "output-pp.txt");
}

void bigSimPP() {
  const int n = 30000;
  std::vector<float> masses = randomMasses(n, {0.002f, 0.002f});
  std::vector<Vec3> state =
      randomInitialState(n, {Vec3(15, 15, 15), Vec3(45, 45, 45)},
                         {Vec3(-0.01f, -0.01f, -0.01f), Vec3(0.01f, 0.01f, 0.01f)});

  int simLength = -1;
  float DT = 1;
  float G = 4.5e-3f;

  ppMethodLeapfrog(state, masses, simLength, DT, G, "output-pp.txt");
}

void smallSimP3M() {
  const int n = 3;
  std::vector<float> masses = {20, 5, 1e2};
  std::vector<Vec3> state = {
      Vec3(30, 30, 15), Vec3(45, 32, 15),  Vec3(30, 10, 15),  // positions
      Vec3(0.1f, 0, 0), Vec3(-0.3f, 0, 0), Vec3::zero()       // velocities
  };

  int simLength = 100;
  float G = 4.5e-3f;

  auto gridPoints = std::make_tuple(64, 64, 32);
  std::array<int, 3> dims = {std::get<2>(gridPoints), std::get<1>(gridPoints),
                             std::get<0>(gridPoints)};
  FFTWAdapter fftAdapter(dims);
  Grid grid(gridPoints, fftAdapter);
  auto effectiveBoxSize = std::make_tuple(60.0f, 60.0f, 30.0f);
  float H = std::get<0>(effectiveBoxSize) / (std::get<0>(gridPoints) / 2);
  float DT = 1;
  auto externalField = [](Vec3 pos) -> Vec3 { return Vec3::zero(); };
  auto externalPotential = [](Vec3 pos) -> float { return 0; };
  float a = 7.5f;

#ifdef CUDA
  PMMethodGPU pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
                 InterpolationScheme::TSC, FiniteDiffScheme::TWO_POINT,
                 GreensFunction::DISCRETE_LAPLACIAN, 0, gridPoints);
#else
  PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
              InterpolationScheme::TSC, FiniteDiffScheme::TWO_POINT, GreensFunction::S1_OPTIMAL, a,
              grid);
#endif

  float re = 0.7f * a;
  float softeningLength = 0.5f;
  P3MMethod p3m(pm, effectiveBoxSize, re, a, H, softeningLength, CloudShape::S1);
  p3m.run(simLength, true /*diagnostics*/, "output-p3m.txt");
}
#include <iostream>
#include "diskSamplerLinear.h"
#include "externalFields.h"
#include "grid.h"
#include "kissFFTAdapter.h"
#include "pmMethod.h"
#include "ppMethod.h"
#include "utils.h"

std::vector<Vec3> linearSpacing(int n, float dx, Vec3 center) {
  std::vector<Vec3> state(2 * n);

  for (int i = 0; i < n; ++i) {
    state[i] = center + Vec3(i * dx, 0, 0);
    state[n + i] = Vec3();
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

  std::vector<Vec3> state = linearSpacing(n, 0.025f /*dx*/, Vec3(30, 30, 30) /*center*/);

  std::vector<float> masses(n, 0);
  masses[0] = 1.0;

  auto externalField = [](Vec3 pos) -> Vec3 { return Vec3(); };
  auto externalPotential = [](Vec3 pos) -> float { return 0; };

  int gridPoints = 64;
  int dims[] = {gridPoints, gridPoints, gridPoints};
  KissFFTAdapter<float> fftAdapter(dims, 3);
  Grid grid(gridPoints, fftAdapter);
  float effectiveBoxSize = 60;
  float H = effectiveBoxSize / (gridPoints / 2);
  float DT = 1;

  PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
              InterpolationScheme::TSC, FiniteDiffScheme::FOUR_POINT, grid);

  pm.run(2, true /*diagnostics*/, true /*field*/);
}

void galaxySimulation() {
  /*
    units:
    [M] = G(solar mass)
    [L] = kpc
    [T] = Myr
    */
  const int n = int(5e4);
  Vec3 galaxyCenter = Vec3(30, 30, 30);
  float rb = 3.0f;
  float mb = 15.0f;
  float rd = 15.0f;
  float md = 60.0f;
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

  int gridPoints = 64;
  int dims[] = {gridPoints, gridPoints, gridPoints};
  KissFFTAdapter<float> fftAdapter(dims, 3);
  Grid grid(gridPoints, fftAdapter);
  float effectiveBoxSize = 60;
  float H = effectiveBoxSize / (gridPoints / 2);
  float DT = 1;

  PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
              InterpolationScheme::TSC, FiniteDiffScheme::TWO_POINT, grid);

  pm.run(200, true /*diagnostics*/);
}

int main() {
  galaxySimulation();
}

#include <iostream>
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

  std::vector<Vec3> state = linearSpacing(n, 0.075f /*dx*/, Vec3(30, 30, 30) /*center*/);

  std::vector<float> masses(n, 0);
  masses[0] = 1.0f;

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
              InterpolationScheme::TSC, FiniteDiffScheme::TWO_POINT, GreensFunction::S1_OPTIMAL,
              grid);

  pm.run(2, true /*diagnostics*/, true /*field*/);
}

void galaxySimulationPM() {
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
              InterpolationScheme::TSC, FiniteDiffScheme::TWO_POINT,
              GreensFunction::DISCRETE_LAPLACIAN, grid);

  pm.run(200, true /*diagnostics*/);
}

void galaxySimulationP3M() {
  /*
    units:
    [M] = G(solar mass)
    [L] = kpc
    [T] = Myr
    */
  const int n = int(1e4);
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
              InterpolationScheme::TSC, FiniteDiffScheme::TWO_POINT, GreensFunction::S1_OPTIMAL,
              grid);

  P3MMethod p3m(pm, effectiveBoxSize, 3.5f, 4.0f, H);

  p3m.run(50, false /*diagnostics*/);
}

void smallSimPP() {
  const int n = 3;
  std::vector<float> masses(n, 1.0f);
  std::vector<Vec3> state = std::vector<Vec3>({
      Vec3(25, 29.5f, 30), Vec3(35, 30.5f, 30), Vec3(30, 35, 30),  // positions
      Vec3(0.1f, 0, 0), Vec3(-0.1f, 0, 0), Vec3(0, -0.11f, 0)      // velocities
  });

  int simLength = 100;
  float DT = 1;
  float G = 4.5e-3f;

  ppMethod(state, masses, simLength, DT, G, "output-pp.txt");
}

void smallSimP3M() {
  const int n = 2;
  std::vector<float> masses(n, 1.0f);
  std::vector<Vec3> state = std::vector<Vec3>({
      Vec3(27, 29.5, 30), Vec3(33, 30.5, 30),  // positions
      Vec3(0.2f, 0, 0), Vec3(-0.2f, 0, 0)      // velocities
  });

  int simLength = 50;
  float G = 4.5e-3f;

  int gridPoints = 64;
  int dims[] = {gridPoints, gridPoints, gridPoints};
  KissFFTAdapter<float> fftAdapter(dims, 3);
  Grid grid(gridPoints, fftAdapter);
  float effectiveBoxSize = 60;
  float H = effectiveBoxSize / (gridPoints / 2);
  float DT = 1;
  auto externalField = [](Vec3 pos) -> Vec3 { return Vec3(); };
  auto externalPotential = [](Vec3 pos) -> float { return 0; };

  PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
              InterpolationScheme::TSC, FiniteDiffScheme::TWO_POINT, GreensFunction::S1_OPTIMAL,
              grid);

  float a = 4.0f;
  float re = 0.7f * a;
  P3MMethod p3m(pm, effectiveBoxSize, re, a, H);
  p3m.run(simLength, true, false, "output-p3m.txt");
  // pm.run(simLength, true, false, "output-pm.txt");
}

int main() {
  // smallSimPP();
  // smallSimP3M();
  // probeField();
  // galaxySimulationP3M();
  galaxySimulationPM();
}

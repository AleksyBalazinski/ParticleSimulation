#include <iostream>
#include "diskSamplerLinear.h"
#include "grid.h"
#include "kissFFTAdapter.h"
#include "pmMethod.h"
#include "ppMethod.h"
#include "utils.h"

int main() {
  /*
  units:
  [M] = G(solar mass)
  [L] = kpc
  [T] = Myr
  */
  const int n = int(5e4);
  Vec3 galaxyCenter = Vec3(30, 30, 30);
  double rb = 3;
  double mb = 15;
  double rd = 15;
  double md = 60;
  double thickness = 0.0;
  double G = 4.5e-3;

  DiskSamplerLinear diskSampler;
  std::vector<Vec3> state = diskSampler.sample(galaxyCenter, rb, mb, rd, md, thickness, G, n);

  std::vector<double> masses(n, md / n);
  auto externalField = [galaxyCenter, rb, mb, G](Vec3 pos) -> Vec3 {
    return externalFieldBulge(pos, galaxyCenter, rb, mb, G);
  };

  int gridPoints = 64;
  int dims[] = {gridPoints, gridPoints, gridPoints};
  KissFFTAdapter<float> fftAdapter(dims, 3);
  Grid grid(gridPoints, fftAdapter);
  double effectiveBoxSize = 60;
  double H = effectiveBoxSize / (gridPoints / 2);
  double DT = 1;

  PMMethod pm(state, masses, externalField, H, DT, G, InterpolationScheme::CIC,
              FiniteDiffScheme::TWO_POINT, grid);

  pm.run(90);
}

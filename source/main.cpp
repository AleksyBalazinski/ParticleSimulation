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
    return externalFieldBulge(pos, galaxyCenter, rb, mb, G);
  };

  int gridPoints = 64;
  int dims[] = {gridPoints, gridPoints, gridPoints};
  KissFFTAdapter<float> fftAdapter(dims, 3);
  Grid grid(gridPoints, fftAdapter);
  float effectiveBoxSize = 60;
  float H = effectiveBoxSize / (gridPoints / 2);
  float DT = 1;

  PMMethod pm(state, masses, effectiveBoxSize, externalField, H, DT, G, InterpolationScheme::CIC,
              FiniteDiffScheme::TWO_POINT, grid);

  pm.run(200, true /*diagnostics*/);
}

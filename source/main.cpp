#include <iostream>
#include "grid.h"
#include "kissFFTAdapter.h"
#include "pmMethod.h"
#include "ppMethod.h"
#include "utils.h"

int main() {
  const int N = 2;
  std::vector<Vec3> state = {Vec3(10, 10, 10), Vec3(25, 15, 20), Vec3(1, 0.5, 0),
                             Vec3(-1, -0.5, 0)};
  std::vector<double> masses = {50, 40};

  int gridPoints = 64;
  int dims[] = {gridPoints, gridPoints, gridPoints};
  KissFFTAdapter<float> fftAdapter(dims, 3);

  Grid grid(gridPoints, fftAdapter);
  PMMethod pm(state, masses, 1, 0.1, 1, InterpolationScheme::CIC, FiniteDiffScheme::TWO_POINT,
              grid);
  pm.run(200);
  // ppMethod(state, masses, /*sim length*/ 200, /*step size*/ 0.1, /*G*/ 1);
}

#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
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

  PMMethod pm(gridPoints, fftAdapter);
  pm.run(state, masses, 15.0, 0.05, 0.5, InterpolationScheme::CIC);
  // ppMethod(state, masses, 10.0);
}
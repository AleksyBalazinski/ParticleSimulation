#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include "pmMethod.h"
#include "ppMethod.h"
#include "utils.h"

int main() {
  const int N = 2;
  std::vector<Vec3> state = {Vec3(10, 10, 10), Vec3(25, 15, 20), Vec3(1, 0.5, 0),
                             Vec3(-1, -0.5, 0)};
  std::vector<double> masses = {50, 40};

  PMMethod pm(state, masses, 32, 1, 0.002);
  pm.run(10.0);
  // ppMethod(state, masses, 10.0);
}
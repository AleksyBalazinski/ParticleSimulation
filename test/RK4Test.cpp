#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <iostream>
#include <numbers>
#include "RK4Stepper.h"

TEST(rk4stepper, harmonicOscillator) {
  typedef std::vector<double> StateType;
  StateType x = {3.0, 0.0};

  auto harmonicOscillator = [](StateType& x, StateType& dxdt, double t) {
    dxdt[0] = x[1];
    dxdt[1] = -x[0];
  };

  RK4Stepper<double> stepper(harmonicOscillator, 2);
  double stopT = 3 * std::numbers::pi;
  double step = 1.0e-4;
  int stepsCnt = (int)(stopT / step);

  for (int i = 0; i < stepsCnt; i++) {
    stepper.doStep(x, step);
  }
  std::cerr << x[0] << ' ' << x[1] << '\n';
  ASSERT_NEAR(x[0], -3.0, 1.0e-3);
  ASSERT_NEAR(x[1], 0.0, 1.0e-3);
}

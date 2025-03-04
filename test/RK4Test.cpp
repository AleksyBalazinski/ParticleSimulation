#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <iostream>
#include <numbers>
#include "RK4Stepper.h"

TEST(rk4stepper, harmonicOscillator) {
  typedef std::vector<float> StateType;
  StateType x = {3.0f, 0.0f};

  auto harmonicOscillator = [](StateType& x, StateType& dxdt, float t) {
    dxdt[0] = x[1];
    dxdt[1] = -x[0];
  };

  RK4Stepper<float> stepper(harmonicOscillator, 2);
  float stopT = 3 * std::numbers::pi_v<float>;
  float step = 1.0e-4f;
  int stepsCnt = (int)(stopT / step);

  for (int i = 0; i < stepsCnt; i++) {
    stepper.doStep(x, step);
  }
  ASSERT_NEAR(x[0], -3.0f, 1.0e-3f);
  ASSERT_NEAR(x[1], 0.0f, 1.0e-3f);
}

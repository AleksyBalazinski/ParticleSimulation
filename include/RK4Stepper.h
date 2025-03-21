#pragma once
#include <array>
#include <functional>
#include "abstractStepper.h"

template <typename T>
class RK4Stepper : public AbstractStepper<T> {
  using StateType = std::vector<T>;
  using SystemFunc = std::function<void(StateType&, StateType&, float)>;

 public:
  RK4Stepper(SystemFunc system, int systemSize, float t0 = 0)
      : system(system),
        dxdt(systemSize),
        k1(systemSize),
        k2(systemSize),
        k3(systemSize),
        k4(systemSize) {
    currentT = t0;
  }

  void doStep(StateType& x, float dt) override {
    system(x, dxdt, currentT);
    Mul(dxdt, dt, k1);

    system(Add(x, Mul(k1, 0.5f, k2), k2), dxdt, currentT + 0.5f * dt);
    Mul(dxdt, dt, k2);

    system(Add(x, Mul(k2, 0.5f, k3), k3), dxdt, currentT + 0.5f * dt);
    Mul(dxdt, dt, k3);

    system(Add(x, k3, k4), dxdt, currentT + dt);
    Mul(dxdt, dt, k4);

    Mul(k2, 2.0f, k2);
    Mul(k3, 2.0f, k3);
    Add(Add(Add(k1, k2, k2), k3, k3), k4, k4);
    Mul(k4, 1.0f / 6, k4);
    Add(x, k4, x);

    currentT += dt;
  }

 private:
  StateType& Mul(StateType& a, float s, StateType& out) {
    for (int i = 0; i < a.size(); i++) {
      out[i] = s * a[i];
    }

    return out;
  }

  StateType& Add(StateType& a, StateType& b, StateType& out) {
    for (int i = 0; i < a.size(); i++) {
      out[i] = a[i] + b[i];
    }

    return out;
  }

  SystemFunc system;
  StateType dxdt;
  float currentT;

  StateType k1;
  StateType k2;
  StateType k3;
  StateType k4;
};
#include "ppMethod.h"

Vec3 acceleration(int i, const std::vector<Vec3>& r, const std::vector<double>& masses, double G) {
  const int n = (int)masses.size();
  Vec3 acc;
  for (int j = 0; j < n; j++) {
    if (j == i)
      continue;

    auto rij = r[i] - r[j];
    acc += -1.0 * G * masses[j] / (rij.getMagnitude() + 0.01) * rij;
  }

  return acc;
}

typedef std::vector<Vec3> StateType;

void nBody(StateType& x, StateType& dxdt, double t, const std::vector<double>& masses, double G) {
  const int n = (int)masses.size();
  for (int i = 0; i < n; i++) {
    dxdt[i] = x[n + i];
    dxdt[n + i] = acceleration(i, x, masses, G);
  }
}

std::string ppMethod(std::vector<Vec3>& state,
                     std::vector<double>& masses,
                     const double simLengthSeconds,
                     const double stepSize,
                     const double G,
                     const char* filepath,
                     const int frameRate) {
  const int n = (int)masses.size();
  auto system = [&masses, G](StateType& x, StateType& dxdt, double t) {
    nBody(x, dxdt, t, masses, G);
  };
  RK4Stepper<Vec3> stepper(system, 2 * n);

  const double frameLength = 1.0 / frameRate;
  StateRecorder stateRecorder(filepath);
  double curFrameAcc = 0;
  for (double t = 0; t <= simLengthSeconds; t += stepSize) {
    if (curFrameAcc <= 0) {
      stateRecorder.record(state.begin(), state.begin() + n);
      curFrameAcc = frameLength;
    }
    stepper.doStep(state, stepSize);
    curFrameAcc -= stepSize;
  }

  stateRecorder.flush();
  return filepath;
}

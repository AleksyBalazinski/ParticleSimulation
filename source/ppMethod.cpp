#include "ppMethod.h"
#include <cmath>

Vec3 acceleration(int i, const std::vector<Vec3>& r, const std::vector<double>& masses, double G) {
  const int n = (int)masses.size();
  Vec3 acc;
  for (int j = 0; j < n; j++) {
    if (j == i)
      continue;

    auto rij = r[i] - r[j];
    acc += -1.0 * G * masses[j] / std::pow((rij.getMagnitude()), 3) * rij;
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

double totalEnergy(StateType& x, const std::vector<double>& masses, double G) {
  const int n = (int)masses.size();
  double potentialEnergy = 0;
  double kineticEnergy = 0;
  for (int i = 0; i < n; i++) {
    kineticEnergy += 0.5 * masses[i] * x[n + i].getMagnitudeSquared();
    for (int j = i + 1; j < n; j++) {
      auto rij = x[i] - x[j];
      potentialEnergy += (-1) * G * masses[i] * masses[j] / rij.getMagnitude();
    }
  }

  return potentialEnergy + kineticEnergy;
}

Vec3 totalMomentum(StateType& x, const std::vector<double>& masses) {
  const int n = (int)masses.size();
  Vec3 momentum;
  for (int i = 0; i < n; i++) {
    momentum += masses[i] * x[n + i];
  }
  return momentum;
}

std::string ppMethod(std::vector<Vec3>& state,
                     std::vector<double>& masses,
                     const double simLengthSeconds,
                     const double stepSize,
                     const double G,
                     const int frameRate,
                     const char* outPath,
                     const char* energyPath,
                     const char* momentumPath) {
  const int n = (int)masses.size();
  auto system = [&masses, G](StateType& x, StateType& dxdt, double t) {
    nBody(x, dxdt, t, masses, G);
  };
  RK4Stepper<Vec3> stepper(system, 2 * n);

  const double frameLength = 1.0 / frameRate;
  StateRecorder stateRecorder(outPath, energyPath, momentumPath);
  double curFrameAcc = 0;
  for (double t = 0; t <= simLengthSeconds; t += stepSize) {
    if (curFrameAcc <= 0) {
      stateRecorder.recordState(state.begin(), state.begin() + n);
      stateRecorder.recordTotalEnergy(totalEnergy(state, masses, G));
      stateRecorder.recordTotalMomentum(totalMomentum(state, masses));
      curFrameAcc = frameLength;
    }
    stepper.doStep(state, stepSize);
    curFrameAcc -= stepSize;
  }

  return stateRecorder.flush();
}

#include "ppMethod.h"
#include <cmath>
#include "simInfo.h"

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

void setIntegerVelocities(const std::vector<Vec3>& state,
                          const std::vector<double>& masses,
                          double G,
                          double h,
                          std::vector<Vec3>& intVs) {
  int n = (int)intVs.size();
  for (size_t i = 0; i < n; i++) {
    intVs[i] = state[n + i] + 0.5 * h * acceleration(i, state, masses, G);
  }
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

  // set v_(1/2)
  // from this point on state[n + i] holds velocities at half-step
  // for (int i = 0; i < n; i++) {
  //   state[n + i] += 0.5 * stepSize * acceleration(i, state, masses, G);
  // }
  // std::vector<Vec3> velocities(n);  // velocities at integer step (needed only for display)

  for (double t = 0; t <= simLengthSeconds; t += stepSize) {
    // for (int i = 0; i < n; i++) {
    //   state[i] += stepSize * state[n + i];
    // }
    if (curFrameAcc <= 0) {
      // setIntegerVelocities(state, masses, G, stepSize, velocities);
      stateRecorder.recordState(state.begin(), state.begin() + n);
      // stateRecorder.recordTotalEnergy(totalEnergy(state.begin(), state.begin() + n,
      //                                             velocities.begin(), velocities.end(), masses,
      //                                             G));
      stateRecorder.recordTotalEnergy(
          totalEnergy(state.begin(), state.begin() + n, state.begin() + n, state.end(), masses, G));
      // stateRecorder.recordTotalMomentum(
      //     totalMomentum(velocities.begin(), velocities.end(), masses));
      stateRecorder.recordTotalMomentum(totalMomentum(state.begin() + n, state.end(), masses));
      curFrameAcc = frameLength;
    }
    stepper.doStep(state, stepSize);
    // for (int i = 0; i < n; i++) {
    //   state[n + i] += stepSize * acceleration(i, state, masses, G);
    // }

    curFrameAcc -= stepSize;
  }

  return stateRecorder.flush();
}

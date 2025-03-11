#include "ppMethod.h"
#include <cmath>
#include "RK4Stepper.h"
#include "simInfo.h"
#include "stateRecorder.h"

// #define LEAPFROG

Vec3 acceleration(int i, const std::vector<Vec3>& r, const std::vector<float>& masses, float G) {
  const int n = (int)masses.size();
  Vec3 acc;
  for (int j = 0; j < n; j++) {
    if (j == i)
      continue;

    auto rij = r[i] - r[j];
    acc += -1.0f * G * masses[j] / std::powf((rij.getMagnitude()), 3) * rij;
  }

  return acc;
}

typedef std::vector<Vec3> StateType;

void nBody(StateType& x, StateType& dxdt, float t, const std::vector<float>& masses, float G) {
  const int n = (int)masses.size();
  for (int i = 0; i < n; i++) {
    dxdt[i] = x[n + i];
    dxdt[n + i] = acceleration(i, x, masses, G);
  }
}

void setIntegerStepVelocities(const std::vector<Vec3>& state,
                              const std::vector<float>& masses,
                              float G,
                              float h,
                              std::vector<Vec3>& intVs) {
  int n = (int)intVs.size();
  for (int i = 0; i < n; i++) {
    intVs[i] = state[n + i] + 0.5f * h * acceleration(i, state, masses, G);
  }
}

std::string ppMethod(std::vector<Vec3>& state,
                     std::vector<float>& masses,
                     const int simLength,
                     const float stepSize,
                     const float G,
                     const char* positionsPath,
                     const char* energyPath,
                     const char* momentumPath) {
  const int n = (int)masses.size();
  StateRecorder stateRecorder(positionsPath, energyPath, momentumPath, "");
#ifdef LEAPFROG
  // set v_(1/2)
  // from this point on state[n + i] holds velocities at half-step
  for (int i = 0; i < n; i++) {
    state[n + i] += 0.5 * stepSize * acceleration(i, state, masses, G);
  }
  std::vector<Vec3> velocities(n);  // velocities at integer step (needed only for display)
#else
  auto system = [&masses, G](StateType& x, StateType& dxdt, float t) {
    nBody(x, dxdt, t, masses, G);
  };

  RK4Stepper<Vec3> stepper(system, 2 * n);
#endif
  for (float t = 0; t <= simLength; ++t) {
#ifdef LEAPFROG
    for (int i = 0; i < n; i++) {
      state[i] += stepSize * state[n + i];
    }
#endif

#ifdef LEAPFROG
    setIntegerStepVelocities(state, masses, G, stepSize, velocities);
    stateRecorder.recordEnergy(potentialEnergy(state.begin(), state.begin() + n, masses, G),
                               kineticEnergy(velocities.begin(), velocities.end(), masses, G));
#endif
    stateRecorder.recordEnergy(
        SimInfo::potentialEnergy(state.begin(), state.begin() + n, masses, G),
        SimInfo::kineticEnergy(state.begin() + n, state.end(), masses, G));
    stateRecorder.recordPositions(state.begin(), state.begin() + n);

#ifdef LEAPFROG
    stateRecorder.recordTotalMomentum(totalMomentum(state.begin() + n, state.end(), masses));
#endif
    stateRecorder.recordTotalMomentum(
        SimInfo::totalMomentum(state.begin() + n, state.end(), masses));

#ifdef LEAPFROG
    for (int i = 0; i < n; i++) {
      state[n + i] += stepSize * acceleration(i, state, masses, G);
    }
#else
    stepper.doStep(state, stepSize);
#endif
  }

  return stateRecorder.flush();
}

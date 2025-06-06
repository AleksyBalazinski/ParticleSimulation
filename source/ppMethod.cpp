#include "ppMethod.h"
#include <cmath>
#include <iostream>
#include "RK4Stepper.h"
#include "leapfrog.h"
#include "simInfo.h"

Vec3 acceleration(int i, const std::vector<Vec3>& r, const std::vector<float>& masses, float G) {
  const int n = (int)masses.size();
  Vec3 acc;
  for (int j = 0; j < n; j++) {
    if (j == i)
      continue;

    Vec3 rij = r[i] - r[j];
    float eps = 0;
    acc += -1.0f * G * masses[j] / std::powf(rij.getMagnitudeSquared() + eps * eps, 1.5f) * rij;
  }

  return acc;
}

void updateAccelerations(std::vector<Particle>& particles, float G, float eps) {
  for (auto& p : particles) {
    p.acceleration = Vec3(0, 0, 0);
  }

  for (int i = 0; i < particles.size(); ++i) {
    for (int j = i + 1; j < particles.size(); ++j) {
      Vec3 rij = particles[i].position - particles[j].position;
      float denom = std::powf(rij.getMagnitudeSquared() + eps * eps, 1.5f);
      particles[i].acceleration += -1 * G * particles[j].mass / denom * rij;
      particles[j].acceleration += G * particles[i].mass / denom * rij;
    }
  }
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

std::string ppMethodRK4(std::vector<Vec3>& state,
                        std::vector<float>& masses,
                        const int simLength,
                        const float stepSize,
                        const float G,
                        StateRecorder& stateRecorder) {
  const int n = (int)masses.size();

  auto system = [&masses, G](StateType& x, StateType& dxdt, float t) {
    nBody(x, dxdt, t, masses, G);
  };

  RK4Stepper<Vec3> stepper(system, 2 * n);

  for (float t = 0; t <= simLength; ++t) {
    stateRecorder.recordEnergy(
        SimInfo::potentialEnergy(state.begin(), state.begin() + n, masses, G),
        SimInfo::kineticEnergy(state.begin() + n, state.end(), masses, G));
    stateRecorder.recordPositions(state.begin(), state.begin() + n);

    stateRecorder.recordTotalMomentum(
        SimInfo::totalMomentum(state.begin() + n, state.end(), masses));

    stepper.doStep(state, stepSize);
  }

  return stateRecorder.flush();
}

std::string ppMethodLeapfrog(const std::vector<Vec3>& state,
                             const std::vector<float>& masses,
                             const int simLength,
                             const float stepSize,
                             const float G,
                             const float softeningLength,
                             StateRecorder& stateRecorder,
                             bool recordField) {
  const int N = int(masses.size());

  std::vector<Particle> particles;
  for (int i = 0; i < N; ++i) {
    particles.emplace_back(state[i], state[N + i], masses[i]);
  }

  updateAccelerations(particles, G, softeningLength);
  setHalfStepVelocities(particles, stepSize);

  for (int t = 0; t <= simLength; ++t) {
    std::cout << "progress: " << float(t) / simLength << '\r';
    std::cout.flush();

    if (recordField) {
      stateRecorder.recordField(particles, 1, 1);
    }

    updatePositions(particles, stepSize);
    stateRecorder.recordPositions(particles);

    setIntegerStepVelocities(particles);
    stateRecorder.recordTotalMomentum(SimInfo::totalMomentum(particles));

    updateAccelerations(particles, G, softeningLength);
    updateVelocities(particles, stepSize);
  }

  return stateRecorder.flush();
}

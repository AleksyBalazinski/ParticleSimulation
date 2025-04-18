#include "unitConversions.h"
#include <algorithm>
#include <execution>
#include <numbers>
#include "particle.h"

void stateToCodeUnits(std::vector<Vec3>& state, float H, float DT) {
  int N = (int)state.size() / 2;
  std::transform(std::execution::par_unseq, state.begin(), state.begin() + N, state.begin(),
                 [H](const Vec3& pos) { return positionToCodeUntits(pos, H); });
  std::transform(std::execution::par_unseq, state.begin() + N, state.end(), state.begin() + N,
                 [H, DT](const Vec3& v) { return velocityToCodeUntits(v, H, DT); });
}

void stateToOriginalUnits(std::vector<Vec3>& state, float H, float DT) {
  int N = (int)state.size() / 2;
  std::transform(std::execution::par_unseq, state.begin(), state.begin() + N, state.begin(),
                 [H](const Vec3& pos) { return positionToOriginalUnits(pos, H); });
  std::transform(std::execution::par_unseq, state.begin() + N, state.end(), state.begin() + N,
                 [H, DT](const Vec3& v) { return velocityToOriginalUnits(v, H, DT); });
}

void stateToCodeUnits(std::vector<Particle>& particles, float H, float DT) {
  std::for_each(std::execution::par_unseq, particles.begin(), particles.end(),
                [H, DT](Particle& p) {
                  p.position = positionToCodeUntits(p.position, H);
                  p.velocity = velocityToCodeUntits(p.velocity, H, DT);
                });
}

void stateToOriginalUnits(std::vector<Particle>& particles, float H, float DT) {
  std::for_each(std::execution::par_unseq, particles.begin(), particles.end(),
                [H, DT](Particle& p) {
                  p.position = positionToOriginalUnits(p.position, H);
                  p.velocity = velocityToOriginalUnits(p.velocity, H, DT);
                });
}

void velocitiesToCodeUnits(std::vector<Vec3>& velocities, float H, float DT) {
  std::transform(velocities.begin(), velocities.end(), velocities.begin(),
                 [H, DT](const Vec3& v) { return velocityToCodeUntits(v, H, DT); });
}

void velocitiesToOriginalUnits(std::vector<Vec3>& velocities, float H, float DT) {
  std::transform(velocities.begin(), velocities.end(), velocities.begin(),
                 [H, DT](const Vec3& v) { return velocityToOriginalUnits(v, H, DT); });
}

void integerStepVelocitiesToOriginalUnits(std::vector<Particle>& particles, float H, float DT) {
  std::for_each(std::execution::par_unseq, particles.begin(), particles.end(),
                [H, DT](Particle& p) {
                  p.integerStepVelocity = velocityToOriginalUnits(p.integerStepVelocity, H, DT);
                });
}

void integerStepVelocitiesToCodeUnits(std::vector<Particle>& particles, float H, float DT) {
  std::for_each(std::execution::par_unseq, particles.begin(), particles.end(),
                [H, DT](Particle& p) {
                  p.integerStepVelocity = velocityToCodeUntits(p.integerStepVelocity, H, DT);
                });
}

void massToCodeUnits(std::vector<Particle>& particles, float H, float DT, float G) {
  std::for_each(std::execution::par_unseq, particles.begin(), particles.end(),
                [H, DT, G](Particle& p) { p.mass = massToCodeUnits(p.mass, H, DT, G); });
}

void massToOriginalUnits(std::vector<Particle>& particles, float H, float DT, float G) {
  std::for_each(std::execution::par_unseq, particles.begin(), particles.end(),
                [H, DT, G](Particle& p) { p.mass = massToOriginalUnits(p.mass, H, DT, G); });
}

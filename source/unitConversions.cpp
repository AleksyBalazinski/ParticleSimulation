#include "unitConversions.h"
#include <algorithm>
#include <execution>
#include <numbers>
#include "particle.h"

Vec3 positionToCodeUntits(const Vec3& pos, float H) {
  return pos / H;
}

Vec3 velocityToCodeUntits(const Vec3& v, float H, float DT) {
  return DT * v / H;
}

Vec3 accelerationToCodeUnits(const Vec3& a, float H, float DT) {
  return DT * DT * a / H;
}

Vec3 accelerationToOriginalUnits(const Vec3& a, float H, float DT) {
  return H * a / (DT * DT);
}

float densityToCodeUnits(float density, float DT, float G) {
  return DT * DT * 4 * std::numbers::pi_v<float> * G * density;
}

float densityToOriginalUnits(float density, float DT, float G) {
  return density / (DT * DT * 4 * std::numbers::pi_v<float> * G);
}

void stateToCodeUnits(std::vector<Vec3>& state, float H, float DT) {
  int N = (int)state.size() / 2;
  std::transform(std::execution::par_unseq, state.begin(), state.begin() + N, state.begin(),
                 [H](const Vec3& pos) { return positionToCodeUntits(pos, H); });
  std::transform(std::execution::par_unseq, state.begin() + N, state.end(), state.begin() + N,
                 [H, DT](const Vec3& v) { return velocityToCodeUntits(v, H, DT); });
}

void stateToCodeUnits(std::vector<Particle>& particles, float H, float DT) {
  std::for_each(std::execution::par_unseq, particles.begin(), particles.end(),
                [H, DT](Particle& p) {
                  p.position = positionToCodeUntits(p.position, H);
                  p.velocity = velocityToCodeUntits(p.velocity, H, DT);
                });
}

void velocitiesToCodeUnits(std::vector<Vec3>& velocities, float H, float DT) {
  std::transform(velocities.begin(), velocities.end(), velocities.begin(),
                 [H, DT](const Vec3& v) { return velocityToCodeUntits(v, H, DT); });
}

Vec3 positionToOriginalUnits(const Vec3& pos, float H) {
  return H * pos;
}

Vec3 velocityToOriginalUnits(const Vec3& v, float H, float DT) {
  return H * v / DT;
}

void stateToOriginalUnits(std::vector<Vec3>& state, float H, float DT) {
  int N = (int)state.size() / 2;
  std::transform(std::execution::par_unseq, state.begin(), state.begin() + N, state.begin(),
                 [H](const Vec3& pos) { return positionToOriginalUnits(pos, H); });
  std::transform(std::execution::par_unseq, state.begin() + N, state.end(), state.begin() + N,
                 [H, DT](const Vec3& v) { return velocityToOriginalUnits(v, H, DT); });
}

void stateToOriginalUnits(std::vector<Particle>& particles, float H, float DT) {
  std::for_each(std::execution::par_unseq, particles.begin(), particles.end(),
                [H, DT](Particle& p) {
                  p.position = positionToOriginalUnits(p.position, H);
                  p.velocity = velocityToOriginalUnits(p.velocity, H, DT);
                });
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

float potentialToOriginalUnits(float potential, float H, float DT) {
  return potential * H * H / (DT * DT);
}

#include "leapfrog.h"
#include <algorithm>
#include <execution>

void setHalfStepVelocities(std::vector<Particle>& particles, float dt) {
  std::for_each(std::execution::par_unseq, particles.begin(), particles.end(),
                [dt](Particle& p) { p.velocity += 0.5f * dt * p.acceleration; });
}

void setIntegerStepVelocities(std::vector<Particle>& particles, float dt) {
  std::for_each(std::execution::par, particles.begin(), particles.end(), [dt](Particle& p) {
    p.integerStepVelocity = p.velocity + 0.5f * dt * p.acceleration;
  });
}

void updateVelocities(std::vector<Particle>& particles, float dt) {
  std::for_each(std::execution::par_unseq, particles.begin(), particles.end(),
                [dt](Particle& p) { p.velocity += dt * p.acceleration; });
}

void updatePositions(std::vector<Particle>& particles, float dt) {
  std::for_each(std::execution::par_unseq, particles.begin(), particles.end(),
                [dt](Particle& p) { p.position += dt * p.velocity; });
}

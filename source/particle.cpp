#include "particle.h"

Particle::Particle(Vec3 position, Vec3 velocity, float mass)
    : position(position),
      velocity(velocity),
      mass(mass),
      integerStepVelocity(velocity),
      HOCNext(-1) {}
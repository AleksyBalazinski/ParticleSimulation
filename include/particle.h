#pragma once
#include "vec3.h"

struct Particle {
  Vec3 position;
  Vec3 velocity;
  Vec3 acceleration;
  float mass;
  Vec3 integerStepVelocity;
  Vec3 shortRangeForce;

  Particle(Vec3 position, Vec3 velocity, float mass);
};
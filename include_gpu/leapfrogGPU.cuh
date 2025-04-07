#pragma once

#include "particle.h"

__global__ void setHalfStepVelocities(Particle* particles, int N, float dt = 1.0f);

__global__ void setIntegerStepVelocities(Particle* particles, int N, float dt = 1.0f);

__global__ void updateVelocities(Particle* particles, int N, float dt = 1.0f);

__global__ void updatePositions(Particle* particles, int N, float dt = 1.0f);

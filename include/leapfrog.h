#pragma once

#include <vector>
#include "particle.h"

void setHalfStepVelocities(std::vector<Particle>& particles, float dt = 1.0f);

void setIntegerStepVelocities(std::vector<Particle>& particles, float dt = 1.0f);

void updateVelocities(std::vector<Particle>& particles, float dt = 1.0f);

void updatePositions(std::vector<Particle>& particles, float dt = 1.0f);
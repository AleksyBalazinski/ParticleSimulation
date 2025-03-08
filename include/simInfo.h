#pragma once

#include <vector>
#include "particle.h"
#include "vec3.h"

float potentialEnergy(std::vector<Vec3>::iterator posBegin,
                      std::vector<Vec3>::iterator posEnd,
                      const std::vector<float>& masses,
                      float G);

float potentialEnergy(const std::vector<Particle>& particles, float G);

float kineticEnergy(std::vector<Vec3>::iterator vBegin,
                    std::vector<Vec3>::iterator vEnd,
                    const std::vector<float>& masses,
                    float G);

float kineticEnergy(const std::vector<Particle>& particles, float G);

float totalEnergy(const std::vector<Vec3>& state, const std::vector<float>& masses, float G);

float totalEnergy(std::vector<Vec3>::iterator posBegin,
                  std::vector<Vec3>::iterator posEnd,
                  std::vector<Vec3>::iterator vBegin,
                  std::vector<Vec3>::iterator vEnd,
                  const std::vector<float>& masses,
                  float G);

Vec3 totalMomentum(const std::vector<Vec3>& state, const std::vector<float>& masses);

Vec3 totalMomentum(std::vector<Vec3>::iterator vBegin,
                   std::vector<Vec3>::iterator vEnd,
                   const std::vector<float>& masses);

Vec3 totalMomentum(const std::vector<Particle>& particles);
#pragma once

#include <vector>
#include "vec3.h"

double potentialEnergy(std::vector<Vec3>::iterator posBegin,
                       std::vector<Vec3>::iterator posEnd,
                       const std::vector<double>& masses,
                       double G);

double kineticEnergy(std::vector<Vec3>::iterator vBegin,
                     std::vector<Vec3>::iterator vEnd,
                     const std::vector<double>& masses,
                     double G);

double totalEnergy(const std::vector<Vec3>& state, const std::vector<double>& masses, double G);

double totalEnergy(std::vector<Vec3>::iterator posBegin,
                   std::vector<Vec3>::iterator posEnd,
                   std::vector<Vec3>::iterator vBegin,
                   std::vector<Vec3>::iterator vEnd,
                   const std::vector<double>& masses,
                   double G);

Vec3 totalMomentum(const std::vector<Vec3>& state, const std::vector<double>& masses);

Vec3 totalMomentum(std::vector<Vec3>::iterator vBegin,
                   std::vector<Vec3>::iterator vEnd,
                   const std::vector<double>& masses);
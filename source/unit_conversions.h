#pragma once
#include <vector>
#include "vec3.h"

Vec3 positionToCodeUntits(const Vec3& pos, double H);

Vec3 velocityToCodeUntits(const Vec3& v, double H, double DT);

double densityToCodeUnits(double density, double DT, double G);

void stateToCodeUnits(std::vector<Vec3>& state, double H, double DT);

void velocitiesToCodeUnits(std::vector<Vec3>& velocities, double H, double DT);

Vec3 positionToOriginalUnits(const Vec3& pos, double H);

Vec3 velocityToOriginalUnits(const Vec3& v, double H, double DT);

void stateToOriginalUnits(std::vector<Vec3>& state, double H, double DT);

void velocitiesToOriginalUnits(std::vector<Vec3>& velocities, double H, double DT);
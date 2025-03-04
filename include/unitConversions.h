#pragma once
#include <vector>
#include "vec3.h"

Vec3 positionToCodeUntits(const Vec3& pos, float H);

Vec3 velocityToCodeUntits(const Vec3& v, float H, float DT);

Vec3 accelerationToCodeUnits(const Vec3& a, float H, float DT);

float densityToCodeUnits(float density, float DT, float G);

void stateToCodeUnits(std::vector<Vec3>& state, float H, float DT);

void velocitiesToCodeUnits(std::vector<Vec3>& velocities, float H, float DT);

Vec3 positionToOriginalUnits(const Vec3& pos, float H);

Vec3 velocityToOriginalUnits(const Vec3& v, float H, float DT);

void stateToOriginalUnits(std::vector<Vec3>& state, float H, float DT);

void velocitiesToOriginalUnits(std::vector<Vec3>& velocities, float H, float DT);
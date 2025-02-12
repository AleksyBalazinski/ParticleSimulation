#pragma once

#include <cmath>
#include <functional>
#include <numbers>
#include <random>
#include <vector>
#include "vec3.h"

std::vector<Vec3> randomInitialState(int particlesCnt,
                                     std::pair<Vec3, Vec3> initalPositionRange,
                                     std::pair<Vec3, Vec3> initialVelocityRange);

std::vector<double> randomMasses(int particlesCnt, std::pair<double, double> massRange);

Vec3 externalFieldBulge(Vec3 pos, Vec3 bulge, double rb, double mb, double G);
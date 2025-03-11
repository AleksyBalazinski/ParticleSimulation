#pragma once

#include <functional>
#include <vector>
#include "vec3.h"

std::vector<Vec3> randomInitialState(int particlesCnt,
                                     std::pair<Vec3, Vec3> initalPositionRange,
                                     std::pair<Vec3, Vec3> initialVelocityRange);

std::vector<float> randomMasses(int particlesCnt, std::pair<float, float> massRange);

float newton(std::function<float(float)> f,
             std::function<float(float)> df,
             float guess,
             int maxIter,
             float tolerance);
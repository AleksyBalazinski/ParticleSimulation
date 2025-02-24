#pragma once

#include "utils.cuh"
#include "vec3.cuh"

__device__ Vec3 positionToCodeUnits(Vec3 pos, float H);

__device__ Vec3 velocityToCodeUnits(Vec3 v, float H, float DT);

__device__ Vec3 positionToOrigUnits(Vec3 pos, float H);

__device__ Vec3 velocityToOrigUnits(Vec3 v, float H, float DT);

__device__ float densityToCodeUnits(float density, float DT, float G);

__device__ Vec3 accelerationToCodeUnits(Vec3 a, float H, float DT);

__global__ void stateToCodeUnits(Vec3* positions, Vec3* velocities, float H, float DT, int n);

__global__ void stateToOrigUnits(Vec3* positions, Vec3* velocities, float H, float DT, int n);
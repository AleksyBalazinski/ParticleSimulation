#include "conversions.cuh"

__device__ Vec3 positionToCodeUnits(Vec3 pos, float H) {
  return pos / H;
}

__device__ Vec3 velocityToCodeUnits(Vec3 v, float H, float DT) {
  return DT * v / H;
}

__device__ Vec3 positionToOrigUnits(Vec3 pos, float H) {
  return H * pos;
}

__device__ Vec3 velocityToOrigUnits(Vec3 v, float H, float DT) {
  return H * v / DT;
}

__device__ float densityToCodeUnits(float density, float DT, float G) {
  return DT * DT * 4 * PI * G * density;
}

__device__ Vec3 accelerationToCodeUnits(Vec3 a, float H, float DT) {
  return DT * DT * a / H;
}

__global__ void stateToCodeUnits(Vec3* positions, Vec3* velocities, float H, float DT, int n) {
  for (int idx = threadIdx.x + blockIdx.x * blockDim.x; idx < n; idx += blockDim.x * gridDim.x) {
    positions[idx] = positionToCodeUnits(positions[idx], H);
    velocities[idx] = velocityToCodeUnits(velocities[idx], H, DT);
  }
}

__global__ void stateToOrigUnits(Vec3* positions, Vec3* velocities, float H, float DT, int n) {
  for (int idx = threadIdx.x + blockIdx.x * blockDim.x; idx < n; idx += blockDim.x * gridDim.x) {
    positions[idx] = positionToOrigUnits(positions[idx], H);
    velocities[idx] = velocityToOrigUnits(velocities[idx], H, DT);
  }
}
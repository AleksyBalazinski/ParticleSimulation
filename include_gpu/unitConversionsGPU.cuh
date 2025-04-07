#pragma once
#include "particle.h"
#include "vec3.h"
#include "vec3GPU.cuh"

#define PI 3.141592f

inline __device__ Vec3 positionToCodeUnits(const Vec3& pos, float H) {
  return div(pos, H);
}
inline __device__ Vec3 positionToOriginalUnits(const Vec3& pos, float H) {
  return mul(H, pos);
}

inline __device__ Vec3 velocityToCodeUnits(const Vec3& v, float H, float DT) {
  return mul(DT, div(v, H));
}
inline __device__ Vec3 velocityToOriginalUnits(const Vec3& v, float H, float DT) {
  return mul(H, div(v, DT));
}

inline __device__ Vec3 accelerationToCodeUnits(const Vec3& a, float H, float DT) {
  return mul(DT * DT, div(a, H));
}

inline __device__ Vec3 accelerationToOriginalUnits(const Vec3& a, float H, float DT) {
  return mul(H, div(a, (DT * DT)));
}

inline __device__ float densityToCodeUnits(float density, float DT, float G) {
  return DT * DT * 4 * PI * G * density;
}
inline __device__ float densityToOriginalUnits(float density, float DT, float G) {
  return density / (DT * DT * 4 * PI * G);
}

inline __device__ float potentialToOriginalUnits(float potential, float H, float DT) {
  return potential * H * H / (DT * DT);
}

inline __device__ float massToCodeUnits(float m, float H, float DT, float G) {
  return DT * DT * 4 * PI * G / (H * H * H) * m;
}
inline __device__ float massToOriginalUnits(float m, float H, float DT, float G) {
  return (H * H * H) / (DT * DT * 4 * PI * G) * m;
}

inline __host__ __device__ float lengthToCodeUnits(float x, float H) {
  return x / H;
}

__global__ void stateToCodeUnits(Particle* particles, int N, float H, float DT);
__global__ void stateToOriginalUnits(Particle* particles, int N, float H, float DT);

__global__ void integerStepVelocitiesToOriginalUnits(Particle* particles, int N, float H, float DT);
__global__ void integerStepVelocitiesToCodeUnits(Particle* particles, int N, float H, float DT);

__global__ void massesToCodeUnits(Particle* particles, int N, float H, float DT, float G);
__global__ void massesToOriginalUnits(Particle* particles, int N, float H, float DT, float G);

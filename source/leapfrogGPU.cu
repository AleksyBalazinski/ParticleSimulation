#include "leapfrogGPU.cuh"
#include "vec3GPU.cuh"

__global__ void setHalfStepVelocities(Particle* particles, int N, float dt) {
  for (int idx = threadIdx.x + blockDim.x * blockIdx.x; idx < N; idx += blockDim.x * gridDim.x) {
    increment(particles[idx].velocity, mul(0.5f * dt, particles[idx].acceleration));
  }
}

__global__ void setIntegerStepVelocities(Particle* particles, int N, float dt) {
  for (int idx = threadIdx.x + blockDim.x * blockIdx.x; idx < N; idx += blockDim.x * gridDim.x) {
    particles[idx].integerStepVelocity =
        add(particles[idx].velocity, mul(0.5f * dt, particles[idx].acceleration));
  }
}

__global__ void updateVelocities(Particle* particles, int N, float dt) {
  for (int idx = threadIdx.x + blockDim.x * blockIdx.x; idx < N; idx += blockDim.x * gridDim.x) {
    increment(particles[idx].velocity, mul(dt, particles[idx].acceleration));
  }
}

__global__ void updatePositions(Particle* particles, int N, float dt) {
  for (int idx = threadIdx.x + blockDim.x * blockIdx.x; idx < N; idx += blockDim.x * gridDim.x) {
    increment(particles[idx].position, mul(dt, particles[idx].velocity));
  }
}

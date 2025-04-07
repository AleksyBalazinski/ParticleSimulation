#include "unitConversionsGPU.cuh"

__global__ void stateToCodeUnits(Particle* particles, int N, float H, float DT) {
  for (int idx = threadIdx.x + blockIdx.x * blockDim.x; idx < N; idx += blockDim.x * gridDim.x) {
    particles[idx].position = positionToCodeUnits(particles[idx].position, H);
    particles[idx].velocity = velocityToCodeUnits(particles[idx].velocity, H, DT);
  }
}

__global__ void stateToOriginalUnits(Particle* particles, int N, float H, float DT) {
  for (int idx = threadIdx.x + blockIdx.x * blockDim.x; idx < N; idx += blockDim.x * gridDim.x) {
    particles[idx].position = positionToOriginalUnits(particles[idx].position, H);
    particles[idx].velocity = velocityToOriginalUnits(particles[idx].velocity, H, DT);
  }
}

__global__ void integerStepVelocitiesToOriginalUnits(Particle* particles,
                                                     int N,
                                                     float H,
                                                     float DT) {
  for (int idx = threadIdx.x + blockIdx.x * blockDim.x; idx < N; idx += blockDim.x * gridDim.x) {
    particles[idx].integerStepVelocity =
        velocityToOriginalUnits(particles[idx].integerStepVelocity, H, DT);
  }
}

__global__ void integerStepVelocitiesToCodeUnits(Particle* particles, int N, float H, float DT) {
  for (int idx = threadIdx.x + blockIdx.x * blockDim.x; idx < N; idx += blockDim.x * gridDim.x) {
    particles[idx].integerStepVelocity =
        velocityToOriginalUnits(particles[idx].integerStepVelocity, H, DT);
  }
}

__global__ void massesToCodeUnits(Particle* particles, int N, float H, float DT, float G) {
  for (int idx = threadIdx.x + blockIdx.x * blockDim.x; idx < N; idx += blockDim.x * gridDim.x) {
    particles[idx].mass = massToCodeUnits(particles[idx].mass, H, DT, G);
  }
}

__global__ void massesToOriginalUnits(Particle* particles, int N, float H, float DT, float G) {
  for (int idx = threadIdx.x + blockIdx.x * blockDim.x; idx < N; idx += blockDim.x * gridDim.x) {
    particles[idx].mass = massToOriginalUnits(particles[idx].mass, H, DT, G);
  }
}

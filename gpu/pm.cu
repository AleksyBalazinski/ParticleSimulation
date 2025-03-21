#include <cufft.h>
#include <algorithm>
#include <chrono>
#include <execution>
#include <iostream>
#include "common.h"
#include "conversions.cuh"
#include "disk_sampler_linear.cuh"
#include "external_fields.cuh"
#include "helper_macros.h"
#include "state_recorder.cuh"
#include "utils.cuh"

#define CELLS_CNT ((NG) * (NG) * (NG))
#define BLOCK_SIZE 1024
#define BLOCK_SIZE_X 16
#define BLOCK_SIZE_Y 8
#define BLOCK_SIZE_Z 8
#define NUM_BLOCKS(n) (((BLOCK_SIZE) + (n) - 1) / (BLOCK_SIZE))
#define GRID_IDX(x, y, z) ((x) + (y) * (NG) + (z) * (NG) * (NG))
#define USE_SMEM

float reassignDensityMs = 0;
float forwardFFTMs = 0;
float findFourierPotentialMs = 0;
float inverseFFTMs = 0;
float scaleAfterInverseMs = 0;
float findFieldInCellsMs = 0;
float updateAccelerationsMs = 0;

cudaEvent_t reassignDensityStart, reassignDensityStop;
cudaEvent_t forwardFFTStart, forwardFFTStop;
cudaEvent_t findFourierPotentialStart, findFourierPotentialStop;
cudaEvent_t inverseFFTStart, inverseFFTStop;
cudaEvent_t scaleAfterInverseStart, scaleAfterInverseStop;
cudaEvent_t findFieldInCellsStart, findFieldInCellsStop;
cudaEvent_t updateAccelerationsStart, updateAccelerationsStop;

dim3 block;
dim3 grid;

float pmTimeMs = 0;
float memcpyTimeMs = 0;
float boundsCheckTimeMs = 0;
float recordStateTimeMs = 0;

__device__ inline int mod(int a, int b) {
  return (a % b + b) % b;
}

__device__ inline int modIndex(int i, int j, int k) {
  return mod(i, NG) + mod(j, NG) * NG + mod(k, NG) * NG * NG;
}

bool isWithinBox(const Vec3 pos, float boxSize) {
  return pos.x >= 0 && pos.x <= boxSize && pos.y >= 0 && pos.y <= boxSize && pos.z >= 0 &&
         pos.z <= boxSize;
}

__global__ void reassignDensity(cufftComplex* gridDensity,
                                Vec3* positions,
                                float* masses,
                                float H,
                                float DT,
                                float G) {
  for (int idx = threadIdx.x + blockDim.x * blockIdx.x; idx < N; idx += blockDim.x * gridDim.x) {
    int particleIdx = idx;

    int x = (int)positions[particleIdx].x;
    int y = (int)positions[particleIdx].y;
    int z = (int)positions[particleIdx].z;

    float vol = H * H * H;
    float d = densityToCodeUnits(masses[particleIdx] / vol, DT, G);

    float dx = positions[particleIdx].x - x;
    float dy = positions[particleIdx].y - y;
    float dz = positions[particleIdx].z - z;
    float tx = 1 - dx;
    float ty = 1 - dy;
    float tz = 1 - dz;

    atomicAdd(reinterpret_cast<float*>(gridDensity) + 2 * GRID_IDX(x, y, z), d * tx * ty * tz);
    atomicAdd(reinterpret_cast<float*>(gridDensity) + 2 * GRID_IDX(x + 1, y, z), d * dx * ty * tz);
    atomicAdd(reinterpret_cast<float*>(gridDensity) + 2 * GRID_IDX(x, y + 1, z), d * tx * dy * tz);
    atomicAdd(reinterpret_cast<float*>(gridDensity) + 2 * GRID_IDX(x, y, z + 1), d * tx * ty * dz);

    atomicAdd(reinterpret_cast<float*>(gridDensity) + 2 * GRID_IDX(x + 1, y + 1, z),
              d * dx * dy * tz);
    atomicAdd(reinterpret_cast<float*>(gridDensity) + 2 * GRID_IDX(x + 1, y, z + 1),
              d * dx * ty * dz);
    atomicAdd(reinterpret_cast<float*>(gridDensity) + 2 * GRID_IDX(x, y + 1, z + 1),
              d * tx * dy * dz);

    atomicAdd(reinterpret_cast<float*>(gridDensity) + 2 * GRID_IDX(x + 1, y + 1, z + 1),
              d * dx * dy * dz);
  }
}

__global__ void findFourierPotential(cufftComplex* potentialFourier,
                                     cufftComplex* gridDensityFourier) {
  for (int idx = threadIdx.x + blockDim.x * blockIdx.x; idx < CELLS_CNT;
       idx += blockDim.x * gridDim.x) {
    int kx = idx % NG;
    int ky = (idx / NG) % NG;
    int kz = idx / (NG * NG);
    if (kx == 0 && ky == 0 && kz == 0) {
      return;
    }
    float sx = sinf(PI * kx / NG);
    float sy = sinf(PI * ky / NG);
    float sz = sinf(PI * kz / NG);
    float green = -0.25f / (sx * sx + sy * sy + sz * sz);

    cufftComplex densityFourier = gridDensityFourier[GRID_IDX(kx, ky, kz)];
    potentialFourier[GRID_IDX(kx, ky, kz)] =
        make_cuComplex(green * densityFourier.x, green * densityFourier.y);
  }
}

__global__ void scaleAfterInverse(cufftComplex* gridPotential) {
  for (int idx = threadIdx.x + blockDim.x * blockIdx.x; idx < CELLS_CNT;
       idx += gridDim.x * blockDim.x) {
    gridPotential[idx].x /= CELLS_CNT;
  }
}

__global__ void findFieldInCells(Vec3* gridField, cufftComplex* gridPotential) {
#ifdef USE_SMEM
  __shared__ float smem_potential[BLOCK_SIZE_X + 2][BLOCK_SIZE_Y + 2][BLOCK_SIZE_Z + 2];
#endif
  int idx = threadIdx.x + blockDim.x * blockIdx.x;
  int idy = threadIdx.y + blockDim.y * blockIdx.y;
  int idz = threadIdx.z + blockDim.z * blockIdx.z;

  if ((idx < NG) && (idy < NG) && (idz < NG)) {
#ifdef USE_SMEM
    int x = threadIdx.x + 1;
    int y = threadIdx.y + 1;
    int z = threadIdx.z + 1;

    smem_potential[x][y][z] = gridPotential[GRID_IDX(idx, idy, idz)].x;

    __syncthreads();

    float xl, yl, zl;
    float xr, yr, zr;

    xl = (x == 1) ? gridPotential[modIndex(idx - 1, idy, idz)].x : smem_potential[x - 1][y][z];
    yl = (y == 1) ? gridPotential[modIndex(idx, idy - 1, idz)].x : smem_potential[x][y - 1][z];
    zl = (z == 1) ? gridPotential[modIndex(idx, idy, idz - 1)].x : smem_potential[x][y][z - 1];
    xr = (x == BLOCK_SIZE_X) ? gridPotential[modIndex(idx + 1, idy, idz)].x
                             : smem_potential[x + 1][y][z];
    yr = (y == BLOCK_SIZE_Y) ? gridPotential[modIndex(idx, idy + 1, idz)].x
                             : smem_potential[x][y + 1][z];
    zr = (z == BLOCK_SIZE_Z) ? gridPotential[modIndex(idx, idy, idz + 1)].x
                             : smem_potential[x][y][z + 1];

    gridField[GRID_IDX(idx, idy, idz)].x = -0.5f * (xr - xl);
    gridField[GRID_IDX(idx, idy, idz)].y = -0.5f * (yr - yl);
    gridField[GRID_IDX(idx, idy, idz)].z = -0.5f * (zr - zl);
#else
    float xl, yl, zl;
    float xr, yr, zr;

    xl = (idx == 0) ? gridPotential[modIndex(idx - 1, idy, idz)].x
                    : gridPotential[GRID_IDX(idx - 1, idy, idz)].x;
    yl = (idy == 0) ? gridPotential[modIndex(idx, idy - 1, idz)].x
                    : gridPotential[GRID_IDX(idx, idy - 1, idz)].x;
    zl = (idz == 0) ? gridPotential[modIndex(idx, idy, idz - 1)].x
                    : gridPotential[GRID_IDX(idx, idy, idz - 1)].x;
    xr = (idx == NG - 1) ? gridPotential[modIndex(idx + 1, idy, idz)].x
                         : gridPotential[GRID_IDX(idx + 1, idy, idz)].x;
    yr = (idy == NG - 1) ? gridPotential[modIndex(idx, idy + 1, idz)].x
                         : gridPotential[GRID_IDX(idx, idy + 1, idz)].x;
    zr = (idz == NG - 1) ? gridPotential[modIndex(idx, idy, idz + 1)].x
                         : gridPotential[GRID_IDX(idx, idy, idz + 1)].x;

    gridField[GRID_IDX(idx, idy, idz)].x = -0.5f * (xr - xl);
    gridField[GRID_IDX(idx, idy, idz)].y = -0.5f * (yr - yl);
    gridField[GRID_IDX(idx, idy, idz)].z = -0.5f * (zr - zl);
#endif
  }
}

__device__ Vec3 interpolateField(Vec3* gridField, Vec3 position) {
  int xi = (int)position.x;
  int yi = (int)position.y;
  int zi = (int)position.z;
  float dx = position.x - xi;
  float dy = position.y - yi;
  float dz = position.z - zi;
  float tx = 1 - dx;
  float ty = 1 - dy;
  float tz = 1 - dz;

  auto field000 = gridField[GRID_IDX(xi, yi, zi)];
  auto field100 = gridField[GRID_IDX(xi + 1, yi, zi)];
  auto field010 = gridField[GRID_IDX(xi, yi + 1, zi)];
  auto field001 = gridField[GRID_IDX(xi, yi, zi + 1)];
  auto field110 = gridField[GRID_IDX(xi + 1, yi + 1, zi)];
  auto field101 = gridField[GRID_IDX(xi + 1, yi, zi + 1)];
  auto field011 = gridField[GRID_IDX(xi, yi + 1, zi + 1)];
  auto field111 = gridField[GRID_IDX(xi + 1, yi + 1, zi + 1)];

  float resX = tx * ty * tz * field000.x + dx * ty * tz * field100.x + tx * dy * tz * field010.x +
               tx * ty * dz * field001.x + dx * dy * tz * field110.x + dx * ty * dz * field101.x +
               tx * dy * dz * field011.x + dx * dy * dz * field111.x;

  float resY = tx * ty * tz * field000.y + dx * ty * tz * field100.y + tx * dy * tz * field010.y +
               tx * ty * dz * field001.y + dx * dy * tz * field110.y + dx * ty * dz * field101.y +
               tx * dy * dz * field011.y + dx * dy * dz * field111.y;

  float resZ = tx * ty * tz * field000.z + dx * ty * tz * field100.z + tx * dy * tz * field010.z +
               tx * ty * dz * field001.z + dx * dy * tz * field110.z + dx * ty * dz * field101.z +
               tx * dy * dz * field011.z + dx * dy * dz * field111.z;

  return Vec3(resX, resY, resZ);
}

__global__ void updateAccelerations(Vec3* accelerations,
                                    Vec3* positions,
                                    Vec3* gridField,
                                    SphRadDecrFieldParams bulgeParams,
                                    float G,
                                    float H,
                                    float DT) {
  for (int idx = threadIdx.x + blockDim.x * blockIdx.x; idx < N; idx += blockDim.x * gridDim.x) {
    Vec3 intField = interpolateField(gridField, positions[idx]);
    Vec3 extField = accelerationToCodeUnits(
        sphRadDecrField(positionToOrigUnits(positions[idx], H), bulgeParams, G), H, DT);

    accelerations[idx].x = intField.x + extField.x;
    accelerations[idx].y = intField.y + extField.y;
    accelerations[idx].z = intField.z + extField.z;
  }
}

__global__ void updateVelocities(Vec3* velocities, Vec3* accelerations) {
  for (int idx = threadIdx.x + blockDim.x * blockIdx.x; idx < N; idx += blockDim.x * gridDim.x) {
    velocities[idx].x += accelerations[idx].x;
    velocities[idx].y += accelerations[idx].y;
    velocities[idx].z += accelerations[idx].z;
  }
}

void pmMethodStep(Vec3* d_accelerations,
                  cufftComplex* d_gridDensity,
                  cufftComplex* d_gridDensityFourier,
                  cufftComplex* d_gridPotential,
                  cufftComplex* d_gridPotentialFourier,
                  Vec3* d_gridField,
                  Vec3* d_positions,
                  float* d_masses,
                  SphRadDecrFieldParams bulgeParams,
                  float H,
                  float DT,
                  float G) {
  cudaMemset(d_gridDensity, 0, CELLS_CNT * sizeof(cufftComplex));
  cudaTime(reassignDensity, reassignDensity<<<NUM_BLOCKS(N), BLOCK_SIZE>>>(
                                d_gridDensity, d_positions, d_masses, H, DT, G));

  cufftHandle plan = 0;
  cufftPlan3d(&plan, NG, NG, NG, CUFFT_C2C);
  cudaTime(forwardFFT, cufftExecC2C(plan, d_gridDensity, d_gridDensityFourier, CUFFT_FORWARD));

  cudaMemset(d_gridPotentialFourier, 0, sizeof(cufftComplex) * CELLS_CNT);
  cudaTime(findFourierPotential, findFourierPotential<<<NUM_BLOCKS(CELLS_CNT), BLOCK_SIZE>>>(
                                     d_gridPotentialFourier, d_gridDensityFourier));

  plan = 0;
  cufftPlan3d(&plan, NG, NG, NG, CUFFT_C2C);
  cudaTime(inverseFFT, cufftExecC2C(plan, d_gridPotentialFourier, d_gridPotential, CUFFT_INVERSE));
  scaleAfterInverse<<<NUM_BLOCKS(CELLS_CNT), BLOCK_SIZE>>>(d_gridPotential);

  cudaTime(findFieldInCells, findFieldInCells<<<grid, block>>>(d_gridField, d_gridPotential));

  cudaTime(updateAccelerations,
           updateAccelerations<<<NUM_BLOCKS(N), BLOCK_SIZE>>>(d_accelerations, d_positions,
                                                              d_gridField, bulgeParams, G, H, DT));

  cudaEventSynchronize(reassignDensityStop);
  cudaEventSynchronize(forwardFFTStop);
  cudaEventSynchronize(findFourierPotentialStop);
  cudaEventSynchronize(inverseFFTStop);
  cudaEventSynchronize(findFieldInCellsStop);
  cudaEventSynchronize(updateAccelerationsStop);

  float milliseconds;
  cudaAccTime(milliseconds, reassignDensity);
  cudaAccTime(milliseconds, forwardFFT);
  cudaAccTime(milliseconds, findFourierPotential);
  cudaAccTime(milliseconds, inverseFFT);
  cudaAccTime(milliseconds, findFieldInCells);
  cudaAccTime(milliseconds, updateAccelerations);
}

__global__ void setHalfVelocities(Vec3* velocities, Vec3* accelerations) {
  for (int idx = threadIdx.x + blockDim.x * blockIdx.x; idx < N; idx += blockDim.x * gridDim.x) {
    velocities[idx] += 0.5f * accelerations[idx];
  }
}

__global__ void updatePositions(Vec3* positions, Vec3* velocities) {
  for (int idx = threadIdx.x + blockDim.x * blockIdx.x; idx < N; idx += blockDim.x * gridDim.x) {
    positions[idx] += velocities[idx];
  }
}

inline float toMs(const std::chrono::nanoseconds delta) {
  return delta.count() * 1e-6f;
}

void pmMethod(std::vector<Vec3>& state,
              const std::vector<float>& masses,
              float effectiveBoxSize,
              float H,
              float DT,
              SphRadDecrFieldParams bulgeParams,
              float G,
              int simLength) {
  block.x = BLOCK_SIZE_X;
  block.y = BLOCK_SIZE_Y;
  block.z = BLOCK_SIZE_Z;
  grid.x = (64 + block.x - 1) / block.x;
  grid.y = (64 + block.y - 1) / block.y;
  grid.z = (64 + block.z - 1) / block.z;

  cudaEventCreate(&reassignDensityStart);
  cudaEventCreate(&reassignDensityStop);

  cudaEventCreate(&forwardFFTStart);
  cudaEventCreate(&forwardFFTStop);

  cudaEventCreate(&findFourierPotentialStart);
  cudaEventCreate(&findFourierPotentialStop);

  cudaEventCreate(&inverseFFTStart);
  cudaEventCreate(&inverseFFTStop);

  cudaEventCreate(&scaleAfterInverseStart);
  cudaEventCreate(&scaleAfterInverseStop);

  cudaEventCreate(&findFieldInCellsStart);
  cudaEventCreate(&findFieldInCellsStop);

  cudaEventCreate(&updateAccelerationsStart);
  cudaEventCreate(&updateAccelerationsStop);

  StateRecorder stateRecorder("output_gpu.txt", "a", "b");

  float* d_masses;
  Vec3* d_positions;
  Vec3* d_velocities;
  Vec3* d_accelerations;

  cufftComplex* d_gridDensity;
  cufftComplex* d_gridDensityFourier;
  cufftComplex* d_gridPotential;
  cufftComplex* d_gridPotentialFourier;
  Vec3* d_gridField;

  cudaMalloc(&d_masses, N * sizeof(float));
  cudaMalloc(&d_positions, N * sizeof(Vec3));
  cudaMalloc(&d_velocities, N * sizeof(Vec3));
  cudaMalloc(&d_accelerations, N * sizeof(Vec3));

  cudaMalloc(&d_gridDensity, CELLS_CNT * sizeof(cufftComplex));
  cudaMalloc(&d_gridDensityFourier, CELLS_CNT * sizeof(cufftComplex));
  cudaMalloc(&d_gridPotential, CELLS_CNT * sizeof(cufftComplex));
  cudaMalloc(&d_gridPotentialFourier, CELLS_CNT * sizeof(cufftComplex));
  cudaMalloc(&d_gridField, CELLS_CNT * sizeof(Vec3));

  cudaMemcpy(d_positions, state.data(), N * sizeof(Vec3), cudaMemcpyHostToDevice);
  cudaMemcpy(d_velocities, state.data() + N, N * sizeof(Vec3), cudaMemcpyHostToDevice);
  cudaMemcpy(d_masses, masses.data(), N * sizeof(float), cudaMemcpyHostToDevice);

  stateToCodeUnits<<<NUM_BLOCKS(N), BLOCK_SIZE>>>(d_positions, d_velocities, H, DT, N);

  pmMethodStep(d_accelerations, d_gridDensity, d_gridDensityFourier, d_gridPotential,
               d_gridPotentialFourier, d_gridField, d_positions, d_masses, bulgeParams, H, DT, G);

  setHalfVelocities<<<NUM_BLOCKS(N), BLOCK_SIZE>>>(d_velocities, d_accelerations);

  for (int t = 0; t <= simLength; ++t) {
    cudaDeviceSynchronize();
    std::cout << "progress: " << float(t) / simLength << '\r';
    std::cout.flush();

    updatePositions<<<NUM_BLOCKS(N), BLOCK_SIZE>>>(d_positions, d_velocities);

    stateToOrigUnits<<<NUM_BLOCKS(N), BLOCK_SIZE>>>(d_positions, d_velocities, H, DT, N);

    hostTime(memcpy,
             cudaMemcpy(state.data(), d_positions, N * sizeof(Vec3), cudaMemcpyDeviceToHost));

    hostTime(recordState, stateRecorder.recordPositions(state.begin(), state.begin() + N));

    auto boundsBegin = std::chrono::steady_clock::now();
    if (std::any_of(
            std::execution::par_unseq, state.begin(), state.begin() + N,
            [effectiveBoxSize](const Vec3& pos) { return !isWithinBox(pos, effectiveBoxSize); })) {
      std::cout << "Particle moved outside the grid.\n";
      break;
    }
    auto boundsEnd = std::chrono::steady_clock::now();
    boundsCheckTimeMs += toMs(boundsEnd - boundsBegin);

    stateToCodeUnits<<<NUM_BLOCKS(N), BLOCK_SIZE>>>(d_positions, d_velocities, H, DT, N);

    hostTime(pm, pmMethodStep(d_accelerations, d_gridDensity, d_gridDensityFourier, d_gridPotential,
                              d_gridPotentialFourier, d_gridField, d_positions, d_masses,
                              bulgeParams, H, DT, G));
    updateVelocities<<<NUM_BLOCKS(N), BLOCK_SIZE>>>(d_velocities, d_accelerations);
  }
  cudaDeviceSynchronize();
  stateRecorder.flush();
  std::cout << "reassign density: " << reassignDensityMs << '\n';
  std::cout << "forward FFT: " << forwardFFTMs << '\n';
  std::cout << "find fourier potential: " << findFourierPotentialMs << '\n';
  std::cout << "inverse FFT: " << inverseFFTMs << '\n';
  std::cout << "find field in cells: " << findFieldInCellsMs << '\n';
  std::cout << "update accelerations: " << updateAccelerationsMs << '\n';
  std::cout << "total PM: " << pmTimeMs << '\n';
  std::cout << "memcpy: " << memcpyTimeMs << '\n';
  std::cout << "bounds check: " << boundsCheckTimeMs << '\n';
  std::cout << "record: " << recordStateTimeMs << '\n';

  cudaFree(d_masses);
  cudaFree(d_positions);
  cudaFree(d_velocities);
  cudaFree(d_accelerations);

  cudaFree(d_gridDensity);
  cudaFree(d_gridDensityFourier);
  cudaFree(d_gridPotential);
  cudaFree(d_gridPotentialFourier);
  cudaFree(d_gridField);

  cudaEventDestroy(reassignDensityStart);
  cudaEventDestroy(reassignDensityStop);

  cudaEventDestroy(forwardFFTStart);
  cudaEventDestroy(forwardFFTStop);

  cudaEventDestroy(findFourierPotentialStart);
  cudaEventDestroy(findFourierPotentialStop);

  cudaEventDestroy(inverseFFTStart);
  cudaEventDestroy(inverseFFTStop);

  cudaEventDestroy(scaleAfterInverseStart);
  cudaEventDestroy(scaleAfterInverseStop);

  cudaEventDestroy(findFieldInCellsStart);
  cudaEventDestroy(findFieldInCellsStop);

  cudaEventDestroy(updateAccelerationsStart);
  cudaEventDestroy(updateAccelerationsStop);
}
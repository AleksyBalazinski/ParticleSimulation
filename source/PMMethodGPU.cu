#include <algorithm>
#include <execution>
#include <iostream>
#include <ranges>
#include "PMMethodGPU.h"
#include "complexGPU.cuh"
#include "greensFunctions.h"
#include "leapfrogGPU.cuh"
#include "measureTime.cuh"
#include "simInfo.h"
#include "unitConversionsGPU.cuh"
#include "vec3GPU.cuh"

#define BLOCK_SIZE 1024
// #define BLOCK_SIZE_X 16
// #define BLOCK_SIZE_Y 8
// #define BLOCK_SIZE_Z 8
#define numBlocks(n) (((BLOCK_SIZE) + (n) - 1) / (BLOCK_SIZE))

declareDeviceTimer(spreadMass);
declareDeviceTimer(fftDensity);
declareDeviceTimer(findFourierPotential);
declareDeviceTimer(invFftPotential);
declareDeviceTimer(findFieldInCells);
declareDeviceTimer(updateAccelerations);

declareHostTimer(memcpy);
declareHostTimer(recordPositions);
declareHostTimer(pmTotal);
declareHostTimer(boundsCheck);

PMMethodGPU::PMMethodGPU(const std::vector<Vec3>& state,
                         const std::vector<float>& masses,
                         const std::tuple<float, float, float> effectiveBoxSize,
                         const std::function<Vec3(Vec3)> externalField,
                         const std::function<float(Vec3)> externalPotential,
                         const float H,
                         const float DT,
                         const float G,
                         const InterpolationScheme is,
                         const FiniteDiffScheme fds,
                         const GreensFunction gFunc,
                         const float particleDiameter,
                         std::tuple<int, int, int> gridPoints)
    : effectiveBoxSizeX(std::get<0>(effectiveBoxSize)),
      effectiveBoxSizeY(std::get<1>(effectiveBoxSize)),
      effectiveBoxSizeZ(std::get<2>(effectiveBoxSize)),
      N(int(masses.size())),
      externalField(externalField),
      externalPotential(externalPotential),
      H(H),
      DT(DT),
      G(G),
      is(is),
      fds(fds),
      gFunc(gFunc),
      particleDiameter(lengthToCodeUnits(particleDiameter, H)),
      grid(gridPoints),
      greensFunction(std::get<0>(gridPoints) * std::get<1>(gridPoints) * std::get<2>(gridPoints)),
      gridDensity(std::get<0>(gridPoints) * std::get<1>(gridPoints) * std::get<2>(gridPoints)),
      gridPotential(std::get<0>(gridPoints) * std::get<1>(gridPoints) * std::get<2>(gridPoints)) {
  for (int i = 0; i < N; ++i) {
    this->particles.emplace_back(state[i], state[N + i], masses[i]);
  }

  cudaMalloc(&d_particles, N * sizeof(Particle));
  cudaMemcpy(d_particles, particles.data(), N * sizeof(Particle), cudaMemcpyHostToDevice);
}

PMMethodGPU::~PMMethodGPU() {
  cudaFree(d_particles);
  grid.freeGrid();
}

std::string PMMethodGPU::run(const int simLength,
                             bool collectDiagnostics,
                             bool recordField,
                             const char* positionsPath,
                             const char* energyPath,
                             const char* momentumPath,
                             const char* expectedMomentumPath,
                             const char* angularMomentumPath,
                             const char* fieldPath) {
  createCudaEvents(spreadMass);
  createCudaEvents(fftDensity);
  createCudaEvents(findFourierPotential);
  createCudaEvents(invFftPotential);
  createCudaEvents(findFieldInCells);
  createCudaEvents(updateAccelerations);

  StateRecorder stateRecorder(positionsPath, N, simLength + 1, energyPath, momentumPath,
                              expectedMomentumPath, angularMomentumPath, fieldPath);
  SimInfo simInfo;

  if (collectDiagnostics) {
    simInfo.setInitialMomentum(particles);
  }

  stateToCodeUnits<<<numBlocks(N), BLOCK_SIZE>>>(d_particles, N, H, DT);
  massesToCodeUnits<<<numBlocks(N), BLOCK_SIZE>>>(d_particles, N, H, DT, G);
  cudaDeviceSynchronize();

  initGreensFunction();
  pmMethodStep();

  setHalfStepVelocities<<<numBlocks(N), BLOCK_SIZE>>>(d_particles, N);

  for (int t = 0; t <= simLength; ++t) {
    cudaDeviceSynchronize();
    std::cout << "progress: " << float(t) / simLength << '\r';
    std::cout.flush();

    updatePositions<<<numBlocks(N), BLOCK_SIZE>>>(d_particles, N);

    if (collectDiagnostics) {
      setIntegerStepVelocities<<<numBlocks(N), BLOCK_SIZE>>>(d_particles, N);
      integerStepVelocitiesToOriginalUnits<<<numBlocks(N), BLOCK_SIZE>>>(d_particles, N, H, DT);
    }

    stateToOriginalUnits<<<numBlocks(N), BLOCK_SIZE>>>(d_particles, N, H, DT);

    if (!collectDiagnostics) {
      measureHostTime(memcpy, cudaMemcpy(particles.data(), d_particles, N * sizeof(Particle),
                                         cudaMemcpyDeviceToHost));
      measureHostTime(recordPositions, stateRecorder.recordPositions(particles));
    } else {
      massesToOriginalUnits<<<numBlocks(N), BLOCK_SIZE>>>(d_particles, N, H, DT, G);
      measureHostTime(memcpy, cudaMemcpy(particles.data(), d_particles, N * sizeof(Particle),
                                         cudaMemcpyDeviceToHost));
      measureHostTime(recordPositions, stateRecorder.recordPositions(particles));

      auto expectedMomentum = simInfo.updateExpectedMomentum(totalExternalForceOrigUnits(), DT);
      stateRecorder.recordExpectedMomentum(expectedMomentum);
      stateRecorder.recordTotalMomentum(SimInfo::totalMomentum(particles));
      stateRecorder.recordTotalAngularMomentum(SimInfo::totalAngularMomentum(particles));
      cudaMemcpy(gridDensity.data(), grid.d_density, grid.getLength() * sizeof(cufftComplex),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(gridPotential.data(), grid.d_potential, grid.getLength() * sizeof(cufftComplex),
                 cudaMemcpyDeviceToHost);
      auto pe = SimInfo::potentialEnergy(gridDensity, gridPotential, particles, externalPotential,
                                         H, DT, G);
      auto ke = SimInfo::kineticEnergy(particles);
      stateRecorder.recordEnergy(pe, ke);
    }
    if (recordField) {
      stateRecorder.recordField(particles, H, DT);
    }

    bool escaped = false;
    measureHostTime(boundsCheck, escaped = escapedComputationalBox());
    if (escaped) {
      std::cout << "Particle moved outside the computational box.\n";
      break;
    }

    stateToCodeUnits<<<numBlocks(N), BLOCK_SIZE>>>(d_particles, N, H, DT);
    if (collectDiagnostics) {
      massesToCodeUnits<<<numBlocks(N), BLOCK_SIZE>>>(d_particles, N, H, DT, G);
      integerStepVelocitiesToCodeUnits<<<numBlocks(N), BLOCK_SIZE>>>(d_particles, N, H, DT);
    }

    measureHostTime(pmTotal, pmMethodStep());

    updateVelocities<<<numBlocks(N), BLOCK_SIZE>>>(d_particles, N);
  }
  cudaDeviceSynchronize();

  printDeviceTime(spreadMass);
  printDeviceTime(fftDensity);
  printDeviceTime(findFourierPotential);
  printDeviceTime(invFftPotential);
  printDeviceTime(findFieldInCells);
  printDeviceTime(updateAccelerations);

  printHostTime(memcpy);
  printHostTime(recordPositions);
  printHostTime(boundsCheck);
  printHostTime(pmTotal);

  destroyCudaEvents(spreadMass);
  destroyCudaEvents(fftDensity);
  destroyCudaEvents(findFourierPotential);
  destroyCudaEvents(invFftPotential);
  destroyCudaEvents(findFieldInCells);
  destroyCudaEvents(updateAccelerations);

  return stateRecorder.flush();
}

void PMMethodGPU::pmMethodStep() {
  measureDeviceTime(spreadMass, spreadMass());
  measureDeviceTime(fftDensity, grid.fftDensity());
  measureDeviceTime(findFourierPotential, findFourierPotential());
  measureDeviceTime(invFftPotential, grid.invFftPotential());
  measureDeviceTime(findFieldInCells, findFieldInCells());
  measureDeviceTime(updateAccelerations, updateAccelerations());

  accDeviceTime(spreadMass);
  accDeviceTime(fftDensity);
  accDeviceTime(findFourierPotential);
  accDeviceTime(invFftPotential);
  accDeviceTime(findFieldInCells);
  accDeviceTime(updateAccelerations);
}

bool PMMethodGPU::escapedComputationalBox() {
  return std::any_of(std::execution::par_unseq, particles.begin(), particles.end(),
                     [this](const Particle& p) { return !isWithingBox(p.position); });
}

Vec3 PMMethodGPU::totalExternalForceOrigUnits() {
  Vec3 f;
  std::for_each(std::execution::seq, particles.begin(), particles.end(),
                [this, &f](const Particle& p) { f += p.mass * externalField(p.position); });

  return f;
}

void PMMethodGPU::initGreensFunction() {
  auto dimsTriple = grid.getGridPoints();
  auto dims = std::make_tuple(dimsTriple.x, dimsTriple.y, dimsTriple.z);
  auto gridRange = std::ranges::views::iota(0, grid.getLength());
  std::for_each(std::execution::par, gridRange.begin(), gridRange.end(), [dims, this](int idx) {
    auto [kx, ky, kz] = grid.indexTripleFromFlat(idx);

    std::complex<float> G;
    if (gFunc == GreensFunction::DISCRETE_LAPLACIAN) {
      G = GreenDiscreteLaplacian(kx, ky, kz, dims);
    } else if (gFunc == GreensFunction::S1_OPTIMAL) {
      if (is == InterpolationScheme::TSC) {
        G = GreenOptimalTSC(kx, ky, kz, dims, particleDiameter, CloudShape::S1, fds);
      } else {
        throw std::invalid_argument("not implemented");
      }
    } else if (gFunc == GreensFunction::S2_OPTIMAL) {
      if (is == InterpolationScheme::TSC) {
        G = GreenOptimalTSC(kx, ky, kz, dims, particleDiameter, CloudShape::S2, fds);
      } else {
        throw std::invalid_argument("not implemented");
      }
    } else {
      throw std::invalid_argument("not implemented");
    }

    greensFunction[idx] = G;
  });

  cudaMemcpy(grid.d_greensFunction, greensFunction.data(), grid.getLength() * sizeof(cufftComplex),
             cudaMemcpyHostToDevice);
}

void PMMethodGPU::copyParticlesDeviceToHost() {
  cudaMemcpy(particles.data(), d_particles, N * sizeof(Particle), cudaMemcpyDeviceToHost);
}

void PMMethodGPU::copyParticlesHostToDevice() {
  cudaMemcpy(d_particles, particles.data(), N * sizeof(Particle), cudaMemcpyHostToDevice);
}

void PMMethodGPU::copyGridPotentialToHost() {
  cudaMemcpy(gridPotential.data(), grid.d_potential, grid.getLength() * sizeof(cufftComplex),
             cudaMemcpyDeviceToHost);
}

void PMMethodGPU::copyGridDensityToHost() {
  cudaMemcpy(gridDensity.data(), grid.d_density, grid.getLength() * sizeof(cufftComplex),
             cudaMemcpyDeviceToHost);
}

__global__ void spreadMassNGPKernel(GridGPU grid, Particle* particles, int N) {
  for (int idx = threadIdx.x + blockDim.x * blockIdx.x; idx < N; idx += blockDim.x * gridDim.x) {
    const Particle& p = particles[idx];
    int x = (int)std::round(p.position.x);
    int y = (int)std::round(p.position.y);
    int z = (int)std::round(p.position.z);

    grid.assignDensity(x, y, z, p.mass);
  }
}

__global__ void spreadMassCICKernel(GridGPU grid, Particle* particles, int N) {
  for (int idx = threadIdx.x + blockDim.x * blockIdx.x; idx < N; idx += blockDim.x * gridDim.x) {
    const Particle& p = particles[idx];
    int x = (int)p.position.x;
    int y = (int)p.position.y;
    int z = (int)p.position.z;

    float d = p.mass;  // divided by vol = H * H * H = 1 (code units)

    float dx = p.position.x - x;
    float dy = p.position.y - y;
    float dz = p.position.z - z;
    float tx = 1 - dx;
    float ty = 1 - dy;
    float tz = 1 - dz;

    grid.assignDensity(x, y, z, d * tx * ty * tz);
    grid.assignDensity(x + 1, y, z, d * dx * ty * tz);
    grid.assignDensity(x, y + 1, z, d * tx * dy * tz);
    grid.assignDensity(x, y, z + 1, d * tx * ty * dz);

    grid.assignDensity(x + 1, y + 1, z, d * dx * dy * tz);
    grid.assignDensity(x + 1, y, z + 1, d * dx * ty * dz);
    grid.assignDensity(x, y + 1, z + 1, d * tx * dy * dz);

    grid.assignDensity(x + 1, y + 1, z + 1, d * dx * dy * dz);
  }
}

__device__ float TSCAssignmentFunc(float x, int t) {
  if (t == 1) {
    return 0.5f * (0.5f + x) * (0.5f + x);
  }
  if (t == 0) {
    return 0.75f - x * x;
  }
  if (t == -1) {
    return 0.5f * (0.5f - x) * (0.5f - x);
  }
  return 0;
}

__global__ void spreadMassTSCKernel(GridGPU grid, Particle* particles, int N) {
  for (int idx = threadIdx.x + blockDim.x * blockIdx.x; idx < N; idx += blockDim.x * gridDim.x) {
    const Particle& p = particles[idx];
    int x = (int)p.position.x;
    int y = (int)p.position.y;
    int z = (int)p.position.z;

    float d = p.mass;

    float dx = p.position.x - x;
    float dy = p.position.y - y;
    float dz = p.position.z - z;

    for (int t1 = -1; t1 <= 1; ++t1) {
      float T1 = (d / 8) * 2 * TSCAssignmentFunc(dx, t1);
      for (int t2 = -1; t2 <= 1; ++t2) {
        float T2 = T1 * 2 * TSCAssignmentFunc(dy, t2);
        for (int t3 = -1; t3 <= 1; ++t3) {
          float T3 = T2 * 2 * TSCAssignmentFunc(dz, t3);
          grid.assignDensity(x + t1, y + t2, z + t3, T3);
        }
      }
    }
  }
}

void PMMethodGPU::spreadMass() {
  grid.clearDensity();
  switch (is) {
    case InterpolationScheme::NGP:
      spreadMassNGPKernel<<<numBlocks(N), BLOCK_SIZE>>>(grid, d_particles, N);
      break;
    case InterpolationScheme::CIC:
      spreadMassCICKernel<<<numBlocks(N), BLOCK_SIZE>>>(grid, d_particles, N);
      break;
    case InterpolationScheme::TSC:
      spreadMassTSCKernel<<<numBlocks(N), BLOCK_SIZE>>>(grid, d_particles, N);
      break;
    default:
      break;
  }
}

__device__ Vec3 interpolateField(Vec3 position, GridGPU grid, InterpolationScheme is) {
  float x = position.x;
  float y = position.y;
  float z = position.z;

  switch (is) {
    case InterpolationScheme::NGP: {
      int xi = (int)std::round(x);
      int yi = (int)std::round(y);
      int zi = (int)std::round(z);
      return grid.getField(xi, yi, zi);
    }

    case InterpolationScheme::CIC: {
      int xi = (int)x;
      int yi = (int)y;
      int zi = (int)z;
      float dx = x - xi;
      float dy = y - yi;
      float dz = z - zi;
      float tx = 1 - dx;
      float ty = 1 - dy;
      float tz = 1 - dz;

      return add(mul(tx * ty * tz, grid.getField(xi, yi, zi)),
                 mul(dx * ty * tz, grid.getField(xi + 1, yi, zi)),
                 mul(tx * dy * tz, grid.getField(xi, yi + 1, zi)),
                 mul(tx * ty * dz, grid.getField(xi, yi, zi + 1)),
                 mul(dx * dy * tz, grid.getField(xi + 1, yi + 1, zi)),
                 mul(dx * ty * dz, grid.getField(xi + 1, yi, zi + 1)),
                 mul(tx * dy * dz, grid.getField(xi, yi + 1, zi + 1)),
                 mul(dx * dy * dz, grid.getField(xi + 1, yi + 1, zi + 1)));
    }

    case InterpolationScheme::TSC: {
      int xi = (int)x;
      int yi = (int)y;
      int zi = (int)z;
      float dx = x - xi;
      float dy = y - yi;
      float dz = z - zi;

      Vec3 interpolatedField;
      for (int t1 = -1; t1 <= 1; ++t1) {
        float T1 = 2 * TSCAssignmentFunc(dx, t1);
        for (int t2 = -1; t2 <= 1; ++t2) {
          float T2 = T1 * 2 * TSCAssignmentFunc(dy, t2);
          for (int t3 = -1; t3 <= 1; ++t3) {
            float T3 = T2 * 2 * TSCAssignmentFunc(dz, t3);
            increment(interpolatedField,
                      mul(0.125f * T3, grid.getField(xi + t1, yi + t2, zi + t3)));
          }
        }
      }

      return interpolatedField;
    }

    default:
      return Vec3(0, 0, 0);
  }
}

__global__ void findFourierPotentialKernel(GridGPU grid) {
  for (int idx = threadIdx.x + blockDim.x * blockIdx.x; idx < grid.getLength();
       idx += blockDim.x * gridDim.x) {
    auto [kx, ky, kz] = grid.indexTripleFromFlat(idx);

    auto densityFourier = grid.getDensityFourier(kx, ky, kz);

    cufftComplex potentialFourier = mul(densityFourier, grid.getGreensFunction(kx, ky, kz));
    grid.setPotentialFourier(kx, ky, kz, potentialFourier);
  }
}

void PMMethodGPU::findFourierPotential() {
  findFourierPotentialKernel<<<numBlocks(grid.getLength()), BLOCK_SIZE>>>(grid);
}

__device__ Vec3 getFieldInCell(int x, int y, int z, FiniteDiffScheme fds, GridGPU grid) {
  float fieldX, fieldY, fieldZ;

  if (fds == FiniteDiffScheme::TWO_POINT) {
    fieldX = -0.5f * (grid.getPotential(x + 1, y, z) - grid.getPotential(x - 1, y, z));
    fieldY = -0.5f * (grid.getPotential(x, y + 1, z) - grid.getPotential(x, y - 1, z));
    fieldZ = -0.5f * (grid.getPotential(x, y, z + 1) - grid.getPotential(x, y, z - 1));
  } else if (fds == FiniteDiffScheme::FOUR_POINT) {
    fieldX = (-1.0f / 12) * (-grid.getPotential(x + 2, y, z) + 8 * grid.getPotential(x + 1, y, z) -
                             8 * grid.getPotential(x - 1, y, z) + grid.getPotential(x - 2, y, z));
    fieldY = (-1.0f / 12) * (-grid.getPotential(x, y + 2, z) + 8 * grid.getPotential(x, y + 1, z) -
                             8 * grid.getPotential(x, y - 1, z) + grid.getPotential(x, y - 2, z));
    fieldZ = (-1.0f / 12) * (-grid.getPotential(x, y, z + 2) + 8 * grid.getPotential(x, y, z + 1) -
                             8 * grid.getPotential(x, y, z - 1) + grid.getPotential(x, y, z - 2));
  } else {
  }

  return Vec3(fieldX, fieldY, fieldZ);
}

__global__ void findFieldInCellsKernel(GridGPU grid, FiniteDiffScheme fds) {
  for (int idx = threadIdx.x + blockDim.x * blockIdx.x; idx < grid.getLength();
       idx += blockDim.x * gridDim.x) {
    auto [x, y, z] = grid.indexTripleFromFlat(idx);
    auto [fieldX, fieldY, fieldZ] = getFieldInCell(x, y, z, fds, grid);

    Vec3 fieldStrength(fieldX, fieldY, fieldZ);
    grid.assignField(x, y, z, fieldStrength);
  }
}

void PMMethodGPU::findFieldInCells() {
  findFieldInCellsKernel<<<numBlocks(grid.getLength()), BLOCK_SIZE>>>(grid, fds);
}

struct UniDecrFieldParams {
  UniDecrFieldParams(Vec3 center, float r, float m, float G) : center(center), r(r), m(m), G(G) {}
  Vec3 center;
  float r;
  float m;
  float G;
};

__device__ Vec3 sphRadDecrField(Vec3 pos, UniDecrFieldParams params) {
  auto bulge = params.center;
  auto mb = params.m;
  auto rb = params.r;

  Vec3 d = sub(pos, bulge);
  float r = norm3df(d.x, d.y, d.z);
  Vec3 dir = div(d, r);
  float g;
  if (r > rb) {
    g = -params.G * mb / (r * r);
  } else {
    g = -(params.G * mb / (rb * rb * rb)) * r * (4 - 3 * r / rb);
  }

  return mul(g, dir);
}

__global__ void updateAccelerationsKernel(Particle* particles,
                                          int N,
                                          GridGPU grid,
                                          InterpolationScheme is,
                                          UniDecrFieldParams params,
                                          float H,
                                          float DT) {
  for (int idx = threadIdx.x + blockDim.x * blockIdx.x; idx < N; idx += blockDim.x * gridDim.x) {
    auto& p = particles[idx];
    p.acceleration =
        add(interpolateField(p.position, grid, is),
            accelerationToCodeUnits(sphRadDecrField(positionToOriginalUnits(p.position, H), params),
                                    H, DT));
  }
}

void PMMethodGPU::updateAccelerations() {
  UniDecrFieldParams params(Vec3(30, 30, 15), 3.0f, 60.0f, 4.5e-3f);  // TODO
  updateAccelerationsKernel<<<numBlocks(N), BLOCK_SIZE>>>(d_particles, N, grid, is, params, H, DT);
}

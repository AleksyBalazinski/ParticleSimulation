#include "gridGPU.cuh"

GridGPU::GridGPU(std::tuple<int, int, int> gridPoints)
    : gridPointsX(std::get<0>(gridPoints)),
      gridPointsY(std::get<1>(gridPoints)),
      gridPointsZ(std::get<2>(gridPoints)),
      length(gridPointsX * gridPointsY * gridPointsZ),
      fftAdapter(gridPoints) {
  cudaMalloc(&d_field, length * sizeof(Vec3));
  cudaMalloc(&d_density, length * sizeof(cufftComplex));
  cudaMalloc(&d_densityFourier, length * sizeof(cufftComplex));
  cudaMalloc(&d_potential, length * sizeof(cufftComplex));
  cudaMalloc(&d_potentialFourier, length * sizeof(cufftComplex));
  cudaMalloc(&d_greensFunction, length * sizeof(cufftComplex));
}

void GridGPU::freeGrid() {
  fftAdapter.free();

  cudaFree(d_field);
  cudaFree(d_density);
  cudaFree(d_densityFourier);
  cudaFree(d_potential);
  cudaFree(d_potentialFourier);
  cudaFree(d_greensFunction);
}

__device__ void GridGPU::assignDensity(int x, int y, int z, float density) {
  atomicAdd((float*)(d_density) + 2 * getIndx(x, y, z), density);
}

void GridGPU::clearDensity() {
  cudaMemset(d_density, 0, length * sizeof(cufftComplex));
}

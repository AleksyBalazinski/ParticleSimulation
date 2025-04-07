#pragma once

#include <cufft.h>
#include <tuple>
#include "cuFFTAdapter.cuh"
#include "triple.cuh"
#include "vec3.h"

class GridGPU {
 public:
  GridGPU(std::tuple<int, int, int> gridPoints);
  void freeGrid();

  __host__ __device__ Triple<int> indexTripleFromFlat(int flatIndex) const {
    int x = flatIndex % gridPointsX;
    int y = (flatIndex / gridPointsX) % gridPointsY;
    int z = flatIndex / (gridPointsX * gridPointsY);
    return makeTriple(x, y, z);
  }

  __device__ void assignDensity(int x, int y, int z, float density);

  void clearDensity();

  __device__ float getDensity(int x, int y, int z) const { return d_density[getIndx(x, y, z)].x; }

  __device__ void assignField(int x, int y, int z, Vec3 fieldVal) {
    d_field[getIndx(x, y, z)] = fieldVal;
  }

  __device__ Vec3 getField(int x, int y, int z) const { return d_field[getIndx(x, y, z)]; }

  int __host__ __device__ getLength() const { return length; }

  Triple<int> getGridPoints() const { return makeTriple(gridPointsX, gridPointsY, gridPointsZ); }

  void fftDensity() { fftAdapter.fft(d_density, d_densityFourier); }

  void invFftPotential() { fftAdapter.ifft(d_potentialFourier, d_potential); }

  __device__ void setPotentialFourier(int i, int j, int k, cufftComplex value) {
    d_potentialFourier[getIndx(i, j, k)] = value;
  }

  __device__ cufftComplex getDensityFourier(int i, int j, int k) const {
    return d_densityFourier[getIndx(i, j, k)];
  }

  __device__ float getPotential(int i, int j, int k) const {
    return d_potential[getWrappedIndx(i, j, k)].x;
  }

  __device__ cufftComplex getGreensFunction(int i, int j, int k) const {
    return d_greensFunction[getIndx(i, j, k)];
  }

  inline __host__ __device__ int getIndx(int i, int j, int k) const {
    return i + j * gridPointsX + k * gridPointsX * gridPointsY;
  }

  cuComplex* d_greensFunction;
  cuComplex* d_density;
  cuComplex* d_potential;

 private:
  inline __device__ int mod(int a, int b) const { return (a % b + b) % b; }
  inline __device__ int getWrappedIndx(int i, int j, int k) const {
    return mod(i, gridPointsX) + mod(j, gridPointsY) * gridPointsX +
           mod(k, gridPointsZ) * gridPointsX * gridPointsY;
  }

  int gridPointsX;
  int gridPointsY;
  int gridPointsZ;
  int length;

  Vec3* d_field;
  cuComplex* d_densityFourier;
  cuComplex* d_potentialFourier;
  CuFFTAdapter fftAdapter;
};
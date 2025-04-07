#pragma once

#include <complex>
#include <mutex>
#include <tuple>
#include <vector>
#include "fftAdapter.h"
#include "vec3.h"

class Grid {
 public:
  Grid(std::tuple<int, int, int> gridPoints, FFTAdapter<float>& fftAdapter);

  std::tuple<int, int, int> indexTripleFromFlat(int flatIndex) const;

  void assignDensity(int x, int y, int z, float density);

  void clearDensity();

  float getDensity(int x, int y, int z) const;

  void assignField(int x, int y, int z, Vec3 fieldVal);

  Vec3 getField(int x, int y, int z) const;

  int getLength() const { return length; }

  std::tuple<int, int, int> getGridPoints() const {
    return std::make_tuple(gridPointsX, gridPointsY, gridPointsZ);
  }

  const std::vector<std::complex<float>>& fftDensity();

  const std::vector<std::complex<float>>& invFftPotential();

  void setPotentialFourier(int i, int j, int k, std::complex<float> value);

  std::complex<float> getDensityFourier(int i, int j, int k) const;

  float getPotential(int i, int j, int k) const;

  std::complex<float> getGreensFunction(int i, int j, int k) const;

  void setGreensFunction(int i, int j, int k, std::complex<float> value);

 private:
  inline int mod(int a, int b) const { return (a % b + b) % b; }
  inline int getWrappedIndx(int i, int j, int k) const {
    return mod(i, gridPointsX) + mod(j, gridPointsY) * gridPointsX +
           mod(k, gridPointsZ) * gridPointsX * gridPointsY;
  }
  inline int getIndx(int i, int j, int k) const {
    return i + j * gridPointsX + k * gridPointsX * gridPointsY;
  }

  int gridPointsX;
  int gridPointsY;
  int gridPointsZ;
  int length;

  std::vector<Vec3> field;
  std::vector<std::complex<float>> density;
  std::vector<std::mutex> densityMutexes;
  std::vector<std::complex<float>> densityFourier;
  std::vector<std::complex<float>> potential;
  std::vector<std::complex<float>> potentialFourier;
  std::vector<std::complex<float>> greensFunction;

  FFTAdapter<float>& fftAdapter;
};

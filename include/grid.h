#pragma once

#include <complex>
#include <vector>
#include "fftAdapter.h"
#include "vec3.h"

class Grid {
 private:
  int gridPoints;
  int length;
  std::vector<Vec3> field;

  long mod(long a, long b) const { return (a % b + b) % b; }

  long getWrappedIndx(int i, int j, int k, int dim) const;

  long getIndx(int i, int j, int k) const;

  std::vector<std::complex<float>> density;
  std::vector<std::complex<float>> densityFourier;
  std::vector<std::complex<float>> potential;
  std::vector<std::complex<float>> potentialFourier;

  FFTAdapter<float>& fftAdapter;

 public:
  Grid(int gridPoints, FFTAdapter<float>& fftAdapter);

  void assignDensity(int x, int y, int z, float density);

  void clearDensity();

  void assignField(int x, int y, int z, Vec3 fieldVal);

  Vec3 getField(int x, int y, int z);

  int getLength() const { return length; }

  int getGridPoints() const { return gridPoints; }

  const std::vector<std::complex<float>>& fftDensity();

  const std::vector<std::complex<float>>& invFftPotential();

  void setPotentialFourier(int i, int j, int k, std::complex<float> value);

  std::complex<float> getDensityFourier(int i, int j, int k) const;

  float getPotential(int i, int j, int k) const;
};

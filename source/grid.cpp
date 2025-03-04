#include "grid.h"
#include <algorithm>
#include <complex>

long Grid::getIndx(int i, int j, int k, int dim) const {
  return mod(i, dim) * dim * dim + mod(j, dim) * dim + mod(k, dim);
}

Grid::Grid(int gridPoints, FFTAdapter<float>& fftAdapter)
    : gridPoints(gridPoints),
      length(gridPoints * gridPoints * gridPoints),
      field(length),
      density(length),
      densityFourier(length),
      potential(length),
      potentialFourier(length),
      fftAdapter(fftAdapter) {
  std::memset(density.data(), 0, length * sizeof(std::complex<float>));
}

void Grid::assignDensity(int x, int y, int z, float d) {
  density[x * gridPoints * gridPoints + y * gridPoints + z] += d;
}

void Grid::clearDensity() {
  std::memset(density.data(), 0, length * sizeof(std::complex<float>));
}

void Grid::assignField(int x, int y, int z, Vec3 fieldVal) {
  field[x * gridPoints * gridPoints + y * gridPoints + z] = fieldVal;
}

Vec3 Grid::getField(int x, int y, int z) {
  return field[x * gridPoints * gridPoints + y * gridPoints + z];
}

const std::vector<std::complex<float>>& Grid::fftDensity() {
  return fftAdapter.fft(density, densityFourier);
}

const std::vector<std::complex<float>>& Grid::invFftPotential() {
  return fftAdapter.ifft(potentialFourier, potential);
}

void Grid::setPotentialFourier(int i, int j, int k, std::complex<float> value) {
  int idx = i * gridPoints * gridPoints + j * gridPoints + k;
  potentialFourier[idx] = value;
}

std::complex<float> Grid::getDensityFourier(int i, int j, int k) const {
  int idx = i * gridPoints * gridPoints + j * gridPoints + k;
  return densityFourier[idx];
}

float Grid::getPotential(int i, int j, int k) const {
  return potential[getIndx(i, j, k, gridPoints)].real();
}

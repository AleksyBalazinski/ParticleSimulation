#include "grid.h"
#include <algorithm>
#include <complex>

long Grid::getWrappedIndx(int i, int j, int k, int dim) const {
  return mod(i, dim) + mod(j, dim) * dim + mod(k, dim) * dim * dim;
}

long Grid::getIndx(int i, int j, int k) const {
  return i + j * gridPoints + k * gridPoints * gridPoints;
}

Grid::Grid(int gridPoints, FFTAdapter<float>& fftAdapter)
    : gridPoints(gridPoints),
      length(gridPoints * gridPoints * gridPoints),
      field(length),
      density(length),
      densityMutexes(length),
      densityFourier(length),
      potential(length),
      potentialFourier(length),
      greensFunction(length),
      fftAdapter(fftAdapter) {
  std::memset(density.data(), 0, length * sizeof(std::complex<float>));
}

std::tuple<int, int, int> Grid::indexTripleFromFlat(int flatIndex) const {
  int x = flatIndex % gridPoints;
  int y = (flatIndex / gridPoints) % gridPoints;
  int z = flatIndex / (gridPoints * gridPoints);
  return std::make_tuple(x, y, z);
}

void Grid::assignDensity(int x, int y, int z, float d) {
  densityMutexes[getIndx(x, y, z)].lock();
  density[getIndx(x, y, z)] += d;
  densityMutexes[getIndx(x, y, z)].unlock();
}

void Grid::clearDensity() {
  std::memset(density.data(), 0, length * sizeof(std::complex<float>));
}

float Grid::getDensity(int x, int y, int z) const {
  return density[getIndx(x, y, z)].real();
}

void Grid::assignField(int x, int y, int z, Vec3 fieldVal) {
  field[getIndx(x, y, z)] = fieldVal;
}

Vec3 Grid::getField(int x, int y, int z) const {
  return field[getIndx(x, y, z)];
}

const std::vector<std::complex<float>>& Grid::fftDensity() {
  return fftAdapter.fft(density, densityFourier);
}

const std::vector<std::complex<float>>& Grid::invFftPotential() {
  return fftAdapter.ifft(potentialFourier, potential);
}

void Grid::setPotentialFourier(int i, int j, int k, std::complex<float> value) {
  potentialFourier[getIndx(i, j, k)] = value;
}

std::complex<float> Grid::getDensityFourier(int i, int j, int k) const {
  return densityFourier[getIndx(i, j, k)];
}

float Grid::getPotential(int i, int j, int k) const {
  return potential[getWrappedIndx(i, j, k, gridPoints)].real();
}

std::complex<float> Grid::getGreensFunction(int i, int j, int k) const {
  return greensFunction[getIndx(i, j, k)];
}

void Grid::setGreensFunction(int i, int j, int k, std::complex<float> value) {
  greensFunction[getIndx(i, j, k)] = value;
}

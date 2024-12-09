#include "grid.h"
#include <algorithm>

long Grid::getIndx(int i, int j, int k, int dim) const {
  return mod(i, dim) * dim * dim + mod(j, dim) * dim + mod(k, dim);
}

Grid::Grid(int gridPoints)
    : gridPoints(gridPoints),
      length(gridPoints * gridPoints * gridPoints),
      field(length),
      density(length),
      densityFourier(length),
      potential(length),
      potentialFourier(length) {
  std::memset(density.data(), 0, length * sizeof(kiss_fft_cpx));

  int dim = gridPoints;
  int dims[] = {dim, dim, dim};
  cfg = kiss_fftnd_alloc(dims, 3, false, nullptr, nullptr);
  cfgInv = kiss_fftnd_alloc(dims, 3, true, nullptr, nullptr);
}

Grid::~Grid() {
  kiss_fft_free(cfg);
  kiss_fft_free(cfgInv);
}

void Grid::assignDensity(int x, int y, int z, double d) {
  density[x * gridPoints * gridPoints + y * gridPoints + z].r += d;
}

void Grid::clearDensity() {
  std::memset(density.data(), 0, length * sizeof(kiss_fft_cpx));
}

void Grid::assignField(int x, int y, int z, Vec3 fieldVal) {
  field[x * gridPoints * gridPoints + y * gridPoints + z] = fieldVal;
}

Vec3 Grid::getField(int x, int y, int z) {
  return field[x * gridPoints * gridPoints + y * gridPoints + z];
}

const std::vector<kiss_fft_cpx>& Grid::fftDensity() {
  kiss_fftnd(cfg, density.data(), densityFourier.data());
  return densityFourier;
}

const std::vector<kiss_fft_cpx>& Grid::invFftPotential() {
  kiss_fftnd(cfgInv, potentialFourier.data(), potential.data());
  for (int i = 0; i < length; i++) {
    potential[i].r /= length;
    potential[i].i /= length;
  }
  return potential;
}

void Grid::setPotentialFourier(int i, int j, int k, kiss_fft_cpx value) {
  int idx = i * gridPoints * gridPoints + j * gridPoints + k;
  potentialFourier[idx].r = value.r;
  potentialFourier[idx].i = value.i;
}

kiss_fft_cpx Grid::getDensityFourier(int i, int j, int k) const {
  int idx = i * gridPoints * gridPoints + j * gridPoints + k;
  return densityFourier[idx];
}

double Grid::getPotential(int i, int j, int k) const {
  return potential[getIndx(i, j, k, gridPoints)].r;
}

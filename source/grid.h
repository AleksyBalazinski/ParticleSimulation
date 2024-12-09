#pragma once

#include <kiss_fftnd.h>
#include <vector>
#include "vec3.h"

class Grid {
 private:
  int gridPoints;
  int length;
  std::vector<Vec3> field;

  long mod(long a, long b) const { return (a % b + b) % b; }

  long getIndx(int i, int j, int k, int dim) const;

  std::vector<kiss_fft_cpx> density;
  std::vector<kiss_fft_cpx> densityFourier;
  std::vector<kiss_fft_cpx> potential;
  std::vector<kiss_fft_cpx> potentialFourier;

  kiss_fftnd_cfg cfg;
  kiss_fftnd_cfg cfgInv;

 public:
  Grid(int gridPoints);

  ~Grid();

  void assignDensity(int x, int y, int z, double density);

  void clearDensity();

  void assignField(int x, int y, int z, Vec3 fieldVal);

  Vec3 getField(int x, int y, int z);

  int getLength() const { return length; }

  int getGridPoints() const { return gridPoints; }

  const std::vector<kiss_fft_cpx>& fftDensity();

  const std::vector<kiss_fft_cpx>& invFftPotential();

  void setPotentialFourier(int i, int j, int k, kiss_fft_cpx value);

  kiss_fft_cpx getDensityFourier(int i, int j, int k) const;

  double getPotential(int i, int j, int k) const;
};

#pragma once

#include <kiss_fftnd.h>
#include <numbers>
#include <string>
#include <vector>
#include "RK4Stepper.h"
#include "grid.h"
#include "simInfo.h"
#include "stateRecorder.h"
#include "vec3.h"

class PMMethod {
 private:
  std::vector<Vec3>& state;
  std::vector<double>& masses;
  std::vector<Vec3> accelerations;

  int gridPoints;
  double H;
  double DT;
  Grid grid;

  // TODO: this should be stored in the Grid class
  kiss_fft_cpx* density;
  kiss_fft_cpx* densityFourier;
  kiss_fft_cpx* potential;
  kiss_fft_cpx* potentialFourier;

  kiss_fftnd_cfg cfg;
  kiss_fftnd_cfg cfgInv;

  Vec3 positionInCodeUntits(const Vec3& pos);
  Vec3 velocityInCodeUnits(const Vec3& v);
  double densityToCodeUnits(double density, double G);
  void stateToCodeUnits(std::vector<Vec3>& state, int n);
  Vec3 positionInOriginalUnits(const Vec3& pos);
  Vec3 velocityInOriginalUnits(const Vec3& v);
  void stateToOriginalUnits(std::vector<Vec3>& state, int n);

  void reassignDensity(const std::vector<double>& masses, double G);
  Vec3 getFieldAtMeshpoint(double x, double y, double z);

  long mod(long a, long b) { return (a % b + b) % b; }

  int getIndx(int i, int j, int k, int dim) {
    return mod(i, dim) * dim * dim + mod(j, dim) * dim + mod(k, dim);
  }

  void updateAccelerations(double G);

 public:
  PMMethod(std::vector<Vec3>& state,
           std::vector<double>& masses,
           int gridPoints,
           double H,
           double DT);
  ~PMMethod();
  std::string run(const double simLengthSeconds = 10.0,
                  const double G = 1,
                  const int frameRate = 30,
                  const char* outPath = "output.txt",
                  const char* energyPath = "energy.txt",
                  const char* momentumPath = "momentum.txt");
};

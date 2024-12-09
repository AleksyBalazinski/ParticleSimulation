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
  Grid grid;

  void reassignDensity(const std::vector<Vec3>& state,
                       const std::vector<double>& masses,
                       double H,
                       double DT,
                       double G);
  Vec3 getFieldAtMeshpoint(double x, double y, double z);

  void updateAccelerations(std::vector<Vec3>& accelerations,
                           const std::vector<Vec3>& state,
                           const std::vector<double>& masses,
                           double H,
                           double DT,
                           double G);

 public:
  PMMethod(int gridPoints);

  std::string run(std::vector<Vec3>& state,
                  std::vector<double>& masses,
                  const double simLengthSeconds = 10.0,
                  const double stepSize = 0.001,
                  const double cellSize = 1,
                  const double G = 1,
                  const int frameRate = 30,
                  const char* outPath = "output.txt",
                  const char* energyPath = "energy.txt",
                  const char* momentumPath = "momentum.txt");
};

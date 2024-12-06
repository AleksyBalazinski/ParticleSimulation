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
  std::vector<Vec3> velocities;  // velocities at integer step (needed only for display)

  int gridPoints;
  int N;
  double H;
  double DT;
  Grid grid;

  Vec3 positionInCodeUntits(const Vec3& pos);
  Vec3 velocityInCodeUnits(const Vec3& v);
  double densityToCodeUnits(double density, double G);
  void stateToCodeUnits();
  void velocitiesToCodeUnits();
  Vec3 positionInOriginalUnits(const Vec3& pos);
  Vec3 velocityInOriginalUnits(const Vec3& v);
  void stateToOriginalUnits();
  void velocitiesToOriginalUnits();

  void reassignDensity(const std::vector<double>& masses, double G);
  Vec3 getFieldAtMeshpoint(double x, double y, double z);

  void updateAccelerations(double G);

 public:
  PMMethod(std::vector<Vec3>& state,
           std::vector<double>& masses,
           int gridPoints,
           double H,
           double DT);  // TODO: store only the grid stuff in the class, there is no reason to force
                        // the user to create a new instance just because they want to run a
                        // different simulation on the same grid

  std::string run(const double simLengthSeconds = 10.0,
                  const double G = 1,
                  const int frameRate = 30,
                  const char* outPath = "output.txt",
                  const char* energyPath = "energy.txt",
                  const char* momentumPath = "momentum.txt");
};

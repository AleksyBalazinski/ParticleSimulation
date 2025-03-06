#pragma once

#include <functional>
#include <vector>
#include "grid.h"
#include "stateRecorder.h"
#include "vec3.h"

enum class InterpolationScheme { NGP, CIC };
enum class FiniteDiffScheme { TWO_POINT, FOUR_POINT };

class PMMethod {
 private:
  Grid& grid;
  float effectiveBoxSize;
  std::vector<Vec3>& state;
  std::vector<float>& masses;
  std::vector<Vec3> accelerations;
  std::vector<Vec3> intStepVelocities;
  int N;
  std::function<Vec3(Vec3)> externalField;
  float H;
  float DT;
  float G;
  InterpolationScheme is;
  FiniteDiffScheme fds;

  void reassignDensity();

  Vec3 interpolateField(float x, float y, float z);

  void setHalfVelocities();

  void updateVelocities();

  void updatePositions();

  void pmMethodStep();

  void findFourierPotential();

  void findFieldInCells();

  void updateAccelerations();

  bool escapedComputationalBox();

 public:
  PMMethod(std::vector<Vec3>& state,
           std::vector<float>& masses,
           float effectiveBoxSize,
           std::function<Vec3(Vec3)> externalField,
           float H,
           float DT,
           float G,
           InterpolationScheme is,
           FiniteDiffScheme fds,
           Grid& grid);

  std::string run(const int simLength,
                  bool collectDiagnostics = false,
                  const char* positionsPath = "output.txt",
                  const char* energyPath = "energy.txt",
                  const char* momentumPath = "momentum.txt");
};

#pragma once

#include <functional>
#include <vector>
#include "grid.h"
#include "particle.h"
#include "stateRecorder.h"
#include "vec3.h"

enum class InterpolationScheme { NGP, CIC };
enum class FiniteDiffScheme { TWO_POINT, FOUR_POINT };

class PMMethod {
 private:
  Grid& grid;
  float effectiveBoxSize;
  std::vector<Particle> particles;
  int N;
  std::function<Vec3(Vec3)> externalField;
  std::function<float(Vec3)> externalPotential;
  float H;
  float DT;
  float G;
  InterpolationScheme is;
  FiniteDiffScheme fds;

  void spreadMass();

  Vec3 interpolateField(Vec3 position);

  void setHalfVelocities();

  void setIntegerStepVelocities();

  void updateVelocities();

  void updatePositions();

  void pmMethodStep();

  void findFourierPotential();

  void findFieldInCells();

  void updateAccelerations();

  bool escapedComputationalBox();

  Vec3 totalExternalForceOrigUnits();

 public:
  PMMethod(const std::vector<Vec3>& state,
           const std::vector<float>& masses,
           const float effectiveBoxSize,
           const std::function<Vec3(Vec3)> externalField,
           const std::function<float(Vec3)> externalPotential,
           const float H,
           const float DT,
           const float G,
           const InterpolationScheme is,
           const FiniteDiffScheme fds,
           Grid& grid);

  std::string run(const int simLength,
                  bool collectDiagnostics = false,
                  const char* positionsPath = "output.txt",
                  const char* energyPath = "energy.txt",
                  const char* momentumPath = "momentum.txt",
                  const char* expectedMomentumPath = "expected_momentum.txt");
};

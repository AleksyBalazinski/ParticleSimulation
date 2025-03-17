#pragma once

#include <functional>
#include <vector>
#include "grid.h"
#include "particle.h"
#include "stateRecorder.h"
#include "vec3.h"

enum class InterpolationScheme { NGP, CIC, TSC };
enum class FiniteDiffScheme { TWO_POINT, FOUR_POINT };
enum class GreensFunction { DISCRETE_LAPLACIAN, S1_OPTIMAL };

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
  GreensFunction gFunc;

  void spreadMass();

  Vec3 interpolateField(Vec3 position);

  void findFourierPotential();

  void findFieldInCells();

  void updateAccelerations();

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
           const GreensFunction gFunc,
           Grid& grid);

  std::string run(const int simLength,
                  bool collectDiagnostics = false,
                  bool recordField = false,
                  const char* positionsPath = "output.txt",
                  const char* energyPath = "energy.txt",
                  const char* momentumPath = "momentum.txt",
                  const char* expectedMomentumPath = "expected_momentum.txt",
                  const char* fieldPath = "field.txt");

  std::vector<Particle>& getParticles();
  float getH() const;
  float getDT() const;
  float getG() const;
  const Grid& getGrid() const;
  std::function<float(Vec3)> getExternalPotential() const;

  void pmMethodStep();

  bool escapedComputationalBox();

  Vec3 totalExternalForceOrigUnits();
};

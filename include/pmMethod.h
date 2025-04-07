#pragma once

#include <functional>
#include <tuple>
#include <vector>
#include "grid.h"
#include "particle.h"
#include "pmConfig.h"
#include "stateRecorder.h"
#include "vec3.h"

class PMMethod {
 public:
  PMMethod(const std::vector<Vec3>& state,
           const std::vector<float>& masses,
           const std::tuple<float, float, float> effectiveBoxSize,
           const std::function<Vec3(Vec3)> externalField,
           const std::function<float(Vec3)> externalPotential,
           const float H,
           const float DT,
           const float G,
           const InterpolationScheme is,
           const FiniteDiffScheme fds,
           const GreensFunction gFunc,
           const float particleDiameter,
           Grid& grid);

  std::string run(const int simLength,
                  bool collectDiagnostics = false,
                  bool recordField = false,
                  const char* positionsPath = "output.dat",
                  const char* energyPath = "energy.txt",
                  const char* momentumPath = "momentum.txt",
                  const char* expectedMomentumPath = "expected_momentum.txt",
                  const char* angularMomentumPath = "angular_momentum.txt",
                  const char* fieldPath = "field.dat");

  std::vector<Particle>& getParticles() { return particles; };
  float getH() const { return H; }
  float getDT() const { return DT; };
  float getG() const { return G; };
  const Grid& getGrid() const { return grid; };
  std::function<float(Vec3)> getExternalPotential() const { return externalPotential; };

  void pmMethodStep();

  bool escapedComputationalBox();

  Vec3 totalExternalForceOrigUnits();

  void initGreensFunction();

 private:
  void spreadMass();
  Vec3 interpolateField(Vec3 position);
  void findFourierPotential();
  void findFieldInCells();
  void updateAccelerations();
  static float TSCAssignmentFunc(float x, int t);

  Grid& grid;
  float effectiveBoxSizeX;
  float effectiveBoxSizeY;
  float effectiveBoxSizeZ;

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

  float particleDiameter;
};

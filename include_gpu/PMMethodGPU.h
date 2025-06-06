#pragma once

#include <functional>
#include <tuple>
#include <vector>
#include "gridGPU.cuh"
#include "particle.h"
#include "pmConfig.h"
#include "stateRecorder.h"
#include "vec3.h"

class PMMethodGPU {
 public:
  PMMethodGPU(const std::vector<Vec3>& state,
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
              std::tuple<int, int, int> gridPoints);
  ~PMMethodGPU();

  std::string run(StateRecorder& stateRecorder,
                  const int simLength,
                  bool collectDiagnostics = false,
                  bool recordField = false);

  std::vector<Particle>& getParticles() { return particles; };
  float getH() const { return H; }
  float getDT() const { return DT; };
  float getG() const { return G; };
  const std::vector<std::complex<float>>& getGridDensity() const { return gridDensity; }
  const std::vector<std::complex<float>>& getGridPotential() const { return gridPotential; }
  std::function<float(Vec3)> getExternalPotential() const { return externalPotential; };

  void pmMethodStep();

  bool escapedComputationalBox();

  Vec3 totalExternalForceOrigUnits();

  void initGreensFunction();

  void copyParticlesDeviceToHost();
  void copyParticlesHostToDevice();

  void copyGridPotentialToHost();
  void copyGridDensityToHost();

 private:
  void spreadMass();
  void findFourierPotential();
  void findFieldInCells();
  void updateAccelerations();

  bool isWithingBox(const Vec3& pos) {
    return pos.x >= 0 && pos.x <= effectiveBoxSizeX && pos.y >= 0 && pos.y <= effectiveBoxSizeY &&
           pos.z >= 0 && pos.z <= effectiveBoxSizeZ;
  }

  float effectiveBoxSizeX;
  float effectiveBoxSizeY;
  float effectiveBoxSizeZ;

  std::vector<Particle> particles;
  Particle* d_particles;
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

  GridGPU grid;
  std::vector<std::complex<float>> greensFunction;
  std::vector<std::complex<float>> gridDensity;
  std::vector<std::complex<float>> gridPotential;
};
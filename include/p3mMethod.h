#pragma once

#ifdef CUDA
#include "PMMethodGPU.h"
#endif
#include "chainingMesh.h"
#include "greensFunctions.h"
#include "pmMethod.h"
#include "simInfo.h"

class P3MMethod {
 public:
  P3MMethod(
#ifdef CUDA
      PMMethodGPU& pmMethod,
#else
      PMMethod& pmMethod,
#endif
      std::tuple<float, float, float> compBoxSize,
      float cutoffRadius,
      float particleDiameter,
      float H,
      float softeningLength,
      CloudShape cloudShape,
      bool useSRForceTable = true);

  void run(const int simLength,
           bool collectDiagnostics = false,
           const char* positionsPath = "output.dat",
           const char* energyPath = "energy.txt",
           const char* momentumPath = "momentum.txt",
           const char* expectedMomentumPath = "expected_momentum.txt",
           const char* angularMomentumPath = "angular_momentum.txt");

 private:
  void calculateShortRangeForces(std::vector<Particle>& particles);

  float referenceForceS1(float r, float a);
  float referenceForceS2(float r, float a);

  Vec3 shortRangeForce(Vec3 rij, float mi, float mj, float a);
  Vec3 shortRangeForceFromTable(Vec3 rij, float mi, float mj, float a);

  void updateSRForces(int i, int j, int q, int qn, int qnLocal, std::vector<Particle>& particles);
  void initSRForceTable();
  void updateSRForcesThreadJob(int tid, int threadsCnt, std::vector<Particle>& particles);

#ifdef CUDA
  PMMethodGPU& pmMethod;
#else
  PMMethod& pmMethod;
#endif
  ChainingMesh chainingMesh;
  SimInfo simInfo;

  float cutoffRadius;
  float particleDiameter;
  int tabulatedValuesCnt;
  float deltaSquared;
  std::vector<float> FTable;
  const float softeningLength;
  bool useSRForceTable;
  const CloudShape cloudShape;
  std::vector<float> threadTimes;
};
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
      bool useSRForceTable = true,
      bool enableYSorting = true);

  void run(StateRecorder& stateRecorder,
           const int simLength,
           bool collectDiagnostics = false,
           bool recordField = false);

 private:
  void calculateShortRangeForces(std::vector<Particle>& particles);

  float referenceForceS1(float r);
  float referenceForceS2(float r);

  Vec3 shortRangeForce(Vec3 rij, float mi, float mj);
  Vec3 shortRangeForceFromTable(Vec3 rij, float mi, float mj);

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
  bool enableYSorting;
  const CloudShape cloudShape;
  std::vector<float> threadTimes;
};
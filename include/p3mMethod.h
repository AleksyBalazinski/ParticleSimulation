#pragma once

#include "chainingMesh.h"
#include "pmMethod.h"
#include "simInfo.h"

class P3MMethod {
 private:
  PMMethod& pmMethod;
  ChainingMesh chainingMesh;
  SimInfo simInfo;
  float cutoffRadius;
  float particleDiameter;
  int tabulatedValuesCnt;
  float deltaSquared;
  std::vector<float> FTable;
  const float softeningLength;
  bool useSRForceTable;

 public:
  P3MMethod(PMMethod& pmMethod,
            float compBoxSize,
            float cutoffRadius,
            float particleDiameter,
            float H,
            float softeningLength,
            bool useSRForceTable = true);

  void run(const int simLength,
           bool collectDiagnostics = false,
           bool recordField = false,
           const char* positionsPath = "output.txt",
           const char* energyPath = "energy.txt",
           const char* momentumPath = "momentum.txt",
           const char* expectedMomentumPath = "expected_momentum.txt",
           const char* fieldPath = "field.txt");

 private:
  void calculateShortRangeForces(std::vector<Particle>& particles);

  void calculateShortRangeForcesPar(std::vector<Particle>& particles);

  float referenceForceS1(float r, float a);

  Vec3 shortRangeForce(Vec3 rij, float mi, float mj, float a);

  Vec3 shortRangeForceFromTable(Vec3 rij, float mi, float mj, float a);

  void updateShortRangeFoces(int i, int j, int q, int qn, std::vector<Particle>& particles);

  void updateSRForcesThreadSafe(int i,
                                int j,
                                int q,
                                int qn,
                                int qnLocal,
                                std::vector<Particle>& particles);

  void initSRForceTable();

  void updateSRForcesThreadJob(int tid, int threadsCnt, std::vector<Particle>& particles);
};
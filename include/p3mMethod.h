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
  float particleRadius;
  int tabulatedValuesCnt;
  float deltaSquared;
  std::vector<float> FTable;

 public:
  P3MMethod(PMMethod& pmMethod,
            float compBoxSize,
            float cutoffRadius,
            float particleRadius,
            float H);

  void calculateShortRangeForces(std::vector<Particle>& particles);

  void run(const int simLength,
           bool collectDiagnostics = false,
           bool recordField = false,
           const char* positionsPath = "output.txt",
           const char* energyPath = "energy.txt",
           const char* momentumPath = "momentum.txt",
           const char* expectedMomentumPath = "expected_momentum.txt",
           const char* fieldPath = "field.txt");

 private:
  float referenceForce(float r, float a);
};
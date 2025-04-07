#pragma once

#include <fstream>
#include <string>
#include <vector>
#include "particle.h"
#include "vec3.h"

class StateRecorder {
 public:
  StateRecorder(const char* positionsPath,
                int particlesCnt,
                int framesCnt,
                const char* energyPath,
                const char* momentumPath,
                const char* expectedMomentumPath,
                const char* angularMomentumPath,
                const char* fieldPath,
                int maxRecords = 500);
  ~StateRecorder();

  void recordPositions(std::vector<Vec3>::iterator begin, std::vector<Vec3>::iterator end);
  void recordPositions(const std::vector<Particle>& particles);

  void recordEnergy(float pe, float ke);
  void recordTotalMomentum(Vec3 momentum);
  void recordExpectedMomentum(Vec3 expectedMomentum);
  void recordTotalAngularMomentum(Vec3 angularMomentum);

  void recordField(const std::vector<Particle>& particles, float H, float DT);

  std::string flush();

 private:
  void saveIfLimitHit(std::ofstream& of, std::string& str, int& counter);
  void saveIfLimitHitBin(std::ofstream& of, std::vector<Vec3>& buf, int& counter);

  std::vector<Vec3> positionsBuf;
  int particlesCnt;
  int framesCnt;
  std::string energyStr;
  std::string momentumStr;
  std::string expectedMomentumStr;
  std::string angularMomentumStr;
  std::vector<Vec3> fieldBuf;

  std::ofstream positionsFile;
  std::ofstream energyFile;
  std::ofstream momentumFile;
  std::ofstream expectedMomentumFile;
  std::ofstream angularMomentumFile;
  std::ofstream fieldFile;

  const char* energyPath;
  const char* momentumPath;
  const char* expectedMomentumPath;
  const char* angularMomentumPath;
  const char* positionsPath;
  const char* fieldPath;

  int energyRecordsCnt = 0;
  int momentumRecordsCnt = 0;
  int expectedMomentumRecordsCnt = 0;
  int angularMomentumRecordsCnt = 0;
  int positionsRecordsCnt = 0;
  int fieldRecordsCnt = 0;

  int maxRecords;

  int vecBufSize = 150;
  std::unique_ptr<char[]> vecBuf;

  int singleBufSize = 50;
  std::unique_ptr<char[]> singleBuf;
};

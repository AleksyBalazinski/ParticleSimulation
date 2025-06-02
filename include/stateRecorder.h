#pragma once

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include "particle.h"
#include "vec3.h"

class StateRecorder {
 public:
  StateRecorder(int particlesCnt,
                int framesCnt,
                const std::filesystem::path& outputDirPath,
                const char* positionsFile = "positions.dat",
                const char* energyFile = "energy.txt",
                const char* momentumFile = "momentum.txt",
                const char* expectedMomentumFile = "expected_momentum.txt",
                const char* angularMomentumFile = "angular_momentum.txt",
                const char* fieldFile = "field.dat",
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

  std::filesystem::path outputDirPath;
  std::filesystem::path energyPath;
  std::filesystem::path momentumPath;
  std::filesystem::path expectedMomentumPath;
  std::filesystem::path angularMomentumPath;
  std::filesystem::path positionsPath;
  std::filesystem::path fieldPath;

  std::ofstream positionsFile;
  std::ofstream energyFile;
  std::ofstream momentumFile;
  std::ofstream expectedMomentumFile;
  std::ofstream angularMomentumFile;
  std::ofstream fieldFile;

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

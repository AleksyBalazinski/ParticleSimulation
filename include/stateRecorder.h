#pragma once

#include <fstream>
#include <string>
#include <vector>
#include "vec3.h"

class StateRecorder {
 private:
  std::string positionsStr;
  std::string energyStr;
  std::string momentumStr;

  std::ofstream positionsFile;
  std::ofstream energyFile;
  std::ofstream momentumFile;

  const char* energyPath;
  const char* momentumPath;
  const char* positionsPath;

  int energyRecordsCnt = 0;
  int momentumRecordsCnt = 0;
  int positionsRecordsCnt = 0;

  int maxRecords;

  int vecBufSize = 150;
  std::unique_ptr<char[]> vecBuf;

  int singleBufSize = 50;
  std::unique_ptr<char[]> singleBuf;

  void saveIfLimitHit(std::ofstream& of, std::string& str, int& counter);

 public:
  StateRecorder(const char* positionsPath,
                const char* energyPath,
                const char* momentumPath,
                int maxRecords = 500);
  ~StateRecorder();
  void recordPositions(std::vector<Vec3>::iterator begin, std::vector<Vec3>::iterator end);
  void recordEnergy(float pe, float ke);
  void recordTotalMomentum(Vec3 momentum);
  std::string flush();
};

#pragma once

#include <string>
#include <vector>
#include "vec3.h"

class StateRecorder {
 private:
  std::string framesStr;
  std::string energyStr;
  std::string momentumStr;

  const char* energyPath;
  const char* momentumPath;
  const char* outPath;

 public:
  StateRecorder(const char* outPath, const char* energyPath, const char* momentumPath);
  void recordPositions(std::vector<Vec3>::iterator begin, std::vector<Vec3>::iterator end);
  void recordEnergy(double pe, double ke);
  void recordTotalMomentum(Vec3 momentum);
  std::string flush();
};

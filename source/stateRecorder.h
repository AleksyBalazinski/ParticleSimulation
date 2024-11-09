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
  void recordState(std::vector<Vec3>::iterator begin, std::vector<Vec3>::iterator end);
  void recordTotalEnergy(double energy);
  void recordTotalMomentum(Vec3 momentum);
  std::string flush();
};

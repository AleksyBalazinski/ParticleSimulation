#pragma once
#include <fstream>
#include <string>
#include <vector>
#include "vec3.h"

class StateRecorder {
 private:
  std::string framesStr;
  std::string filepath;

 public:
  StateRecorder(std::string filepath);
  void record(std::vector<Vec3>::iterator begin, std::vector<Vec3>::iterator end);
  void flush();
};

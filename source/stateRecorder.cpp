#include "stateRecorder.h"

StateRecorder::StateRecorder(std::string filepath) : filepath(std::move(filepath)) {}

void StateRecorder::record(std::vector<Vec3>::iterator begin, std::vector<Vec3>::iterator end) {
  for (auto it = begin; it != end; ++it) {
    framesStr += it->toString() + '\n';
  }
  framesStr += "\n\n";
}

void StateRecorder::flush() {
  std::ofstream recordFile(filepath, std::ofstream::trunc);
  recordFile << framesStr;
  recordFile.close();
}
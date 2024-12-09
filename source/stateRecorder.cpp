#include "stateRecorder.h"
#include <filesystem>
#include <fstream>

StateRecorder::StateRecorder(const char* outPath, const char* energyPath, const char* momentumPath)
    : outPath(outPath), energyPath(energyPath), momentumPath(momentumPath) {}

void StateRecorder::recordPositions(std::vector<Vec3>::iterator begin,
                                    std::vector<Vec3>::iterator end) {
  for (auto it = begin; it != end; ++it) {
    framesStr += it->toString() + '\n';
  }
  framesStr += "\n\n";
}

void StateRecorder::recordEnergy(double pe, double ke) {
  energyStr += std::to_string(pe) + ' ' + std::to_string(ke);
  energyStr += '\n';
}

void StateRecorder::recordTotalMomentum(Vec3 momentum) {
  momentumStr += momentum.toString();
  momentumStr += '\n';
}

std::string StateRecorder::flush() {
  std::ofstream recordFile(outPath, std::ofstream::trunc);
  std::ofstream energyFile(energyPath, std::ofstream::trunc);
  std::ofstream momentumFile(momentumPath, std::ofstream::trunc);

  recordFile << framesStr;
  energyFile << energyStr;
  momentumFile << momentumStr;

  recordFile.close();
  energyFile.close();
  momentumFile.close();

  std::filesystem::path cwd = std::filesystem::current_path();
  return cwd.string();
}
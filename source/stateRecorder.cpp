#include "stateRecorder.h"
#include <filesystem>

void StateRecorder::saveIfLimitHit(std::ofstream& of, std::string& str, int& counter) {
  if (counter < maxRecords)
    return;

  of << str;
  str.clear();
  counter = 0;
}

StateRecorder::StateRecorder(const char* positionsPath,
                             const char* energyPath,
                             const char* momentumPath,
                             int maxRecords)
    : positionsPath(positionsPath),
      energyPath(energyPath),
      momentumPath(momentumPath),
      maxRecords(maxRecords),
      positionsFile(positionsPath, std::ofstream::trunc),
      energyFile(energyPath, std::ofstream::trunc),
      momentumFile(momentumPath, std::ofstream::trunc) {}

StateRecorder::~StateRecorder() {
  positionsFile.close();
  energyFile.close();
  momentumFile.close();
}

void StateRecorder::recordPositions(std::vector<Vec3>::iterator begin,
                                    std::vector<Vec3>::iterator end) {
  for (auto it = begin; it != end; ++it) {
    positionsStr += it->toString() + '\n';
  }
  positionsStr += "\n\n";
  ++positionsRecordsCnt;
  saveIfLimitHit(positionsFile, positionsStr, positionsRecordsCnt);
}

void StateRecorder::recordEnergy(double pe, double ke) {
  energyStr += std::to_string(pe) + ' ' + std::to_string(ke);
  energyStr += '\n';
  ++energyRecordsCnt;
  saveIfLimitHit(energyFile, energyStr, energyRecordsCnt);
}

void StateRecorder::recordTotalMomentum(Vec3 momentum) {
  momentumStr += momentum.toString();
  momentumStr += '\n';
  ++momentumRecordsCnt;
  saveIfLimitHit(momentumFile, momentumStr, momentumRecordsCnt);
}

std::string StateRecorder::flush() {
  positionsFile << positionsStr;
  energyFile << energyStr;
  momentumFile << momentumStr;

  std::filesystem::path cwd = std::filesystem::current_path();
  return cwd.string();
}
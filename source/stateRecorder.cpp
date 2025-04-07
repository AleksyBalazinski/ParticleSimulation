#include "stateRecorder.h"
#include <filesystem>
#include "unitConversions.h"

StateRecorder::StateRecorder(const char* positionsPath,
                             int particlesCnt,
                             int framesCnt,
                             const char* energyPath,
                             const char* momentumPath,
                             const char* expectedMomentumPath,
                             const char* angularMomentumPath,
                             const char* fieldPath,
                             int maxRecords)
    : positionsPath(positionsPath),
      particlesCnt(particlesCnt),
      framesCnt(framesCnt),
      energyPath(energyPath),
      momentumPath(momentumPath),
      expectedMomentumPath(expectedMomentumPath),
      angularMomentumPath(angularMomentumPath),
      fieldPath(fieldPath),
      maxRecords(maxRecords),
      energyFile(energyPath, std::ofstream::trunc),
      momentumFile(momentumPath, std::ofstream::trunc),
      expectedMomentumFile(expectedMomentumPath, std::ofstream::trunc),
      angularMomentumFile(angularMomentumPath, std::ofstream::trunc),
      vecBuf(new char[vecBufSize]),
      singleBuf(new char[singleBufSize]) {
  {
    std::ofstream clearFile(positionsPath, std::ios::trunc);
  }

  positionsFile.open(positionsPath, std::ios::binary | std::ios::app);
  positionsFile.write(reinterpret_cast<char*>(&this->particlesCnt), sizeof(int));
  positionsFile.write(reinterpret_cast<char*>(&this->framesCnt), sizeof(int));

  {
    std::ofstream clearFile(fieldPath, std::ios::trunc);
  }

  fieldFile.open(fieldPath, std::ios::binary | std::ios::app);
  fieldFile.write(reinterpret_cast<char*>(&this->particlesCnt), sizeof(int));
  fieldFile.write(reinterpret_cast<char*>(&this->framesCnt), sizeof(int));
}

StateRecorder::~StateRecorder() {
  positionsFile.close();
  energyFile.close();
  momentumFile.close();
}

void StateRecorder::recordPositions(std::vector<Vec3>::iterator begin,
                                    std::vector<Vec3>::iterator end) {
  for (auto it = begin; it != end; ++it) {
    positionsBuf.emplace_back(it->x, it->y, it->z);
    ++positionsRecordsCnt;
    saveIfLimitHitBin(positionsFile, positionsBuf, positionsRecordsCnt);
  }
}

void StateRecorder::recordPositions(const std::vector<Particle>& particles) {
  for (const auto& p : particles) {
    positionsBuf.emplace_back(p.position.x, p.position.y, p.position.z);
    ++positionsRecordsCnt;
    saveIfLimitHitBin(positionsFile, positionsBuf, positionsRecordsCnt);
  }
}

void StateRecorder::recordEnergy(float pe, float ke) {
  energyStr += std::to_string(pe) + ' ' + std::to_string(ke);
  energyStr += '\n';
  ++energyRecordsCnt;
  saveIfLimitHit(energyFile, energyStr, energyRecordsCnt);
}

void StateRecorder::recordTotalMomentum(Vec3 momentum) {
  momentumStr += momentum.toString(singleBuf.get(), singleBufSize, vecBuf.get(), vecBufSize);
  momentumStr += '\n';
  ++momentumRecordsCnt;
  saveIfLimitHit(momentumFile, momentumStr, momentumRecordsCnt);
}

void StateRecorder::recordExpectedMomentum(Vec3 expectedMomentum) {
  expectedMomentumStr +=
      expectedMomentum.toString(singleBuf.get(), singleBufSize, vecBuf.get(), vecBufSize);
  expectedMomentumStr += '\n';
  ++expectedMomentumRecordsCnt;
  saveIfLimitHit(expectedMomentumFile, expectedMomentumStr, expectedMomentumRecordsCnt);
}

void StateRecorder::recordTotalAngularMomentum(Vec3 angularMomentum) {
  angularMomentumStr +=
      angularMomentum.toString(singleBuf.get(), singleBufSize, vecBuf.get(), vecBufSize);
  angularMomentumStr += '\n';
  ++angularMomentumRecordsCnt;
  saveIfLimitHit(angularMomentumFile, angularMomentumStr, angularMomentumRecordsCnt);
}

void StateRecorder::recordField(const std::vector<Particle>& particles, float H, float DT) {
  for (const auto& p : particles) {
    fieldBuf.push_back(accelerationToOriginalUnits(p.acceleration, H, DT));
    ++fieldRecordsCnt;
    saveIfLimitHitBin(fieldFile, fieldBuf, fieldRecordsCnt);
  }
}

std::string StateRecorder::flush() {
  positionsFile.write(reinterpret_cast<char*>(positionsBuf.data()),
                      positionsBuf.size() * sizeof(Vec3));
  energyFile << energyStr;
  momentumFile << momentumStr;
  expectedMomentumFile << expectedMomentumStr;
  angularMomentumFile << angularMomentumStr;
  fieldFile.write(reinterpret_cast<char*>(fieldBuf.data()), fieldBuf.size() * sizeof(Vec3));

  std::filesystem::path cwd = std::filesystem::current_path();
  return cwd.string();
}

void StateRecorder::saveIfLimitHit(std::ofstream& of, std::string& str, int& counter) {
  if (counter < maxRecords)
    return;

  of << str;
  str.clear();
  counter = 0;
}

void StateRecorder::saveIfLimitHitBin(std::ofstream& of, std::vector<Vec3>& buf, int& counter) {
  if (counter < maxRecords)
    return;

  of.write(reinterpret_cast<char*>(buf.data()), buf.size() * sizeof(Vec3));
  buf.clear();
  counter = 0;
}

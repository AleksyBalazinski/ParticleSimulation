#pragma once

#include <array>
#include <memory>
#include <tuple>
#include <vector>
#include "particle.h"

class ChainingMesh {
 public:
  ChainingMesh(std::tuple<float, float, float> compBoxSize, float cutoffRadius, float H, int N);

  void fillWithYSorting(std::vector<Particle>& particles);

  std::array<int, 14> getNeighborsAndSelf(int cellIdx) const;

  inline int getParticlesInCell(int cellIdx) { return hoc[cellIdx]; }
  inline int getSize() const { return size; };
  inline std::tuple<int, int, int> getLength() const { return std::make_tuple(Mx, My, Mz); }
  inline int listEnd() const { return -1; }

 private:
  int tripleToFlatIndex(int x, int y, int z) const;
  inline std::tuple<int, int, int> flatToTripleIndex(int idx) const {
    return std::make_tuple(idx % Mx, (idx / Mx) % My, idx / (Mx * My));
  };

  int Mx;
  int My;
  int Mz;
  float HCx;
  float HCy;
  float HCz;
  int size;
  std::vector<int> hoc;
};
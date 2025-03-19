#pragma once

#include <array>
#include <memory>
#include <tuple>
#include <vector>
#include "particle.h"

class ChainingMesh {
 public:
  struct LLNode;

 private:
  int M;
  float HC;
  int size;
  std::vector<LLNode*> hoc;
  std::unique_ptr<LLNode[]> nodePool;

 public:
  struct LLNode {
    int particleId;
    LLNode* next;
    LLNode(int particleId, LLNode* next);
    LLNode() = default;
  };

  ChainingMesh(float compBoxSize, float cutoffRadius, float H, int N);

  void fillWithYSorting(const std::vector<Particle>& particles);

  std::array<int, 14> getNeighborsAndSelf(int cellIdx) const;

  LLNode* getParticlesInCell(int cellIdx);

  int getSize() const;

  int getLength() const;

 private:
  int tripleToFlatIndex(int x, int y, int z) const;

  std::tuple<int, int, int> flatToTripleIndex(int idx) const;
};
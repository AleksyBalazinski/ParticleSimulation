#pragma once

#include <array>
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

 public:
  struct LLNode {
    int particleId;
    LLNode* next;
    LLNode(int particleId, LLNode* next);
  };

  ChainingMesh(float compBoxSize, float cutoffRadius, float H);

  ~ChainingMesh();

  void fill(const std::vector<Particle>& particles);

  void clear();

  std::array<int, 14> getNeighborsAndSelf(int cellIdx);

  LLNode* getParticlesInCell(int cellIdx);

  int getSize();

  int getLength();

 private:
  int tripleToFlatIndex(int x, int y, int z);

  std::tuple<int, int, int> flatToTripleIndex(int idx);
};
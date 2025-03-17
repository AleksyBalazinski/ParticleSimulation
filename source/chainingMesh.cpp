#include "chainingMesh.h"
#include "unitConversions.h"

ChainingMesh::ChainingMesh(float compBoxSize, float cutoffRadius, float H)
    : M(int(compBoxSize / cutoffRadius)),
      HC(lengthToCodeUnits(compBoxSize / M, H)),
      size(M * M * M),
      hoc(size, nullptr) {}

ChainingMesh::~ChainingMesh() {
  clear();
}

void ChainingMesh::fill(const std::vector<Particle>& particles) {
  std::memset(hoc.data(), 0, size * sizeof(LLNode*));

  for (int i = 0; i < particles.size(); ++i) {
    const auto& p = particles[i];
    int cellX = int(p.position.x / HC);
    int cellY = int(p.position.y / HC);
    int cellZ = int(p.position.z / HC);

    int cellIdx = tripleToFlatIndex(cellX, cellY, cellZ);
    LLNode* head = hoc[cellIdx];
    hoc[cellIdx] = new LLNode(i, head);
  }
}

void ChainingMesh::clear() {
  for (int i = 0; i < size; ++i) {
    for (LLNode* node = hoc[i]; node != nullptr;) {
      LLNode* next = node->next;
      delete node;
      node = next;
    }
  }
}

std::array<int, 14> ChainingMesh::getNeighborsAndSelf(int cellIdx) {
  auto [cellX, cellY, cellZ] = flatToTripleIndex(cellIdx);
  std::array<int, 14> neighbors;

  int i = 0;
  for (int t = -1; t <= 1; ++t) {
    for (int s = -1; s <= 1; ++s) {
      neighbors[i++] = tripleToFlatIndex(cellX + t, cellY - 1, cellZ + s);
    }
  }
  for (int s = -1; s <= 1; ++s) {
    neighbors[i++] = tripleToFlatIndex(cellX + s, cellY, cellZ - 1);
  }
  neighbors[12] = tripleToFlatIndex(cellX - 1, cellY, cellZ);
  neighbors[13] = cellIdx;

  return neighbors;
}

ChainingMesh::LLNode* ChainingMesh::getParticlesInCell(int cellIdx) {
  return hoc[cellIdx];
}

int ChainingMesh::getSize() {
  return size;
}

int ChainingMesh::getLength() {
  return M;
}

int ChainingMesh::tripleToFlatIndex(int x, int y, int z) {
  if (x < 0 || y < 0 || z < 0 || x >= M || y >= M || z >= M) {
    return -1;
  }
  return x + y * M + z * M * M;
}

std::tuple<int, int, int> ChainingMesh::flatToTripleIndex(int idx) {
  return std::make_tuple(idx % M, (idx / M) % M, idx / (M * M));
}

ChainingMesh::LLNode::LLNode(int particleId, LLNode* next) : particleId(particleId), next(next) {}

#include "chainingMesh.h"
#include "unitConversions.h"

ChainingMesh::ChainingMesh(float compBoxSize, float cutoffRadius, float H, int N)
    : M(int(compBoxSize / cutoffRadius)),
      HC(lengthToCodeUnits(compBoxSize / M, H)),
      size(M * M * M),
      hoc(size, nullptr),
      nodePool(new LLNode[N]) {}

void ChainingMesh::fillWithYSorting(const std::vector<Particle>& particles) {
  std::memset(hoc.data(), 0, size * sizeof(LLNode*));

  for (int i = 0; i < particles.size(); ++i) {
    const auto& p = particles[i];
    int cellX = int(p.position.x / HC);
    int cellY = int(p.position.y / HC);
    int cellZ = int(p.position.z / HC);

    int cellIdx = tripleToFlatIndex(cellX, cellY, cellZ);
    LLNode* head = hoc[cellIdx];
    if (head == nullptr || particles[head->particleId].position.y > p.position.y) {
      hoc[cellIdx] = new (nodePool.get() + i) LLNode(i, head);
      continue;
    }

    for (LLNode* node = head; node != nullptr; node = node->next) {
      if (node->next == nullptr || particles[node->next->particleId].position.y > p.position.y) {
        node->next = new (nodePool.get() + i) LLNode(i, node->next);
        break;
      }
    }
  }
}

std::array<int, 14> ChainingMesh::getNeighborsAndSelf(int cellIdx) const {
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

int ChainingMesh::getSize() const {
  return size;
}

int ChainingMesh::getLength() const {
  return M;
}

int ChainingMesh::tripleToFlatIndex(int x, int y, int z) const {
  if (x < 0 || y < 0 || z < 0 || x >= M || y >= M || z >= M) {
    return -1;
  }
  return x + y * M + z * M * M;
}

std::tuple<int, int, int> ChainingMesh::flatToTripleIndex(int idx) const {
  return std::make_tuple(idx % M, (idx / M) % M, idx / (M * M));
}

ChainingMesh::LLNode::LLNode(int particleId, LLNode* next) : particleId(particleId), next(next) {}

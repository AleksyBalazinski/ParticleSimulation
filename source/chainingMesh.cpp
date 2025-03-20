#include "chainingMesh.h"
#include "unitConversions.h"

ChainingMesh::ChainingMesh(std::tuple<float, float, float> compBoxSize,
                           float cutoffRadius,
                           float H,
                           int N)
    : Mx(int(std::get<0>(compBoxSize) / cutoffRadius)),
      My(int(std::get<1>(compBoxSize) / cutoffRadius)),
      Mz(int(std::get<2>(compBoxSize) / cutoffRadius)),
      HCx(lengthToCodeUnits(std::get<0>(compBoxSize) / Mx, H)),
      HCy(lengthToCodeUnits(std::get<1>(compBoxSize) / My, H)),
      HCz(lengthToCodeUnits(std::get<2>(compBoxSize) / Mz, H)),
      size(Mx * My * Mz),
      hoc(size, nullptr),
      nodePool(new LLNode[N]) {}

void ChainingMesh::fillWithYSorting(const std::vector<Particle>& particles) {
  std::memset(hoc.data(), 0, size * sizeof(LLNode*));

  for (int i = 0; i < particles.size(); ++i) {
    const auto& p = particles[i];
    int cellX = int(p.position.x / HCx);
    int cellY = int(p.position.y / HCy);
    int cellZ = int(p.position.z / HCz);

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

std::tuple<int, int, int> ChainingMesh::getLength() const {
  return std::make_tuple(Mx, My, Mz);
}

int ChainingMesh::tripleToFlatIndex(int x, int y, int z) const {
  if (x < 0 || y < 0 || z < 0 || x >= Mx || y >= My || z >= Mz) {
    return -1;
  }
  return x + y * Mx + z * Mx * My;
}

std::tuple<int, int, int> ChainingMesh::flatToTripleIndex(int idx) const {
  return std::make_tuple(idx % Mx, (idx / Mx) % My, idx / (Mx * My));
}

ChainingMesh::LLNode::LLNode(int particleId, LLNode* next) : particleId(particleId), next(next) {}

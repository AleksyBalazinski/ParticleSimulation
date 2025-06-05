#include "barnesHut.h"
#include <algorithm>
#include <bitset>
#include <chrono>
#include <cmath>
#include <execution>
#include <iostream>
#include <numeric>
#include <ranges>
#include "simInfo.h"
#include "stateRecorder.h"

namespace BH {

int nodeID = 0;

Node::Node(Vec3 low, float H, Node* parent, Particle* p)
    : low(low), H(H), p(p), COM(), M(0.0f), children({nullptr}), parent(parent) {
  id = nodeID++;
}

Node::~Node() {}

bool nodeContains(const Particle& p, const Node& n) {
  return p.position.x > n.low.x && p.position.x < n.low.x + n.H && p.position.y > n.low.y &&
         p.position.y < n.low.y + n.H && p.position.z > n.low.z && p.position.z < n.low.z + n.H;
}

Tree::Tree(const std::vector<Particle>& particles, Vec3 low, float H, bool useZOrdering)
    : compBoxLow(low), compBoxSize(H), zCodes(particles.size()), lastInserted(nullptr) {
  if (useZOrdering) {
    initZCodes(particles);
    std::vector<int> indices(particles.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(std::execution::par, indices.begin(), indices.end(),
              [this](int i, int j) { return zCodes[i] < zCodes[j]; });

    root = new Node(low, H, nullptr);
    for (int idx : indices) {  // iterate particles in Z-order
      const auto& p = particles[idx];
      if (lastInserted == nullptr) {
        insert(root, p);
      } else {
        auto node = lastInserted;
        while (!nodeContains(p, *node)) {
          node = node->parent;
        }
        insert(node, p);
      }
    }
  } else {
    root = new Node(low, H, nullptr);

    for (const auto& p : particles) {
      insert(root, p);
    }
  }
}

Tree::~Tree() {
  DFSFree(root);
  delete root;
}

void Tree::print() {
  DFSPrint(root);
}

const Node* Tree::getRoot() const {
  return root;
}

uint32_t expandBits3D(uint32_t v) {
  v &= 0x000003ff;
  v = (v | (v << 16)) & 0x030000FF;
  v = (v | (v << 8)) & 0x0300F00F;
  v = (v | (v << 4)) & 0x030C30C3;
  v = (v | (v << 2)) & 0x09249249;
  return v;
}

uint32_t zCode3D(float x, float y, float z, float compBoxSize, Vec3 compBoxLow) {
  uint32_t resolution = 1024;
  x = std::clamp(x - compBoxLow.x, 0.0f, compBoxSize);
  y = std::clamp(y - compBoxLow.y, 0.0f, compBoxSize);
  z = std::clamp(z - compBoxLow.z, 0.0f, compBoxSize);

  uint32_t ix = static_cast<uint32_t>((x / compBoxSize) * (resolution - 1));
  uint32_t iy = static_cast<uint32_t>((y / compBoxSize) * (resolution - 1));
  uint32_t iz = static_cast<uint32_t>((z / compBoxSize) * (resolution - 1));

  return (expandBits3D(iz) << 2) | (expandBits3D(iy) << 1) | expandBits3D(ix);
}

void Tree::initZCodes(const std::vector<Particle>& particles) {
  for (int i = 0; i < particles.size(); ++i) {
    const auto& pos = particles[i].position;
    zCodes[i] = zCode3D(pos.x, pos.y, pos.z, compBoxSize, compBoxLow);
  }
}

Vec3 operator*(const Mat3x3& m, const Vec3& v) {
  return {m(0, 0) * v.x + m(0, 1) * v.y + m(0, 2) * v.z,
          m(1, 0) * v.x + m(1, 1) * v.y + m(1, 2) * v.z,
          m(2, 0) * v.x + m(2, 1) * v.y + m(2, 2) * v.z};
}

Mat3x3 outerProd(const Vec3& u, const Vec3& v) {
  Mat3x3 result;
  result(0, 0) = u.x * v.x;
  result(0, 1) = u.x * v.y;
  result(0, 2) = u.x * v.z;

  result(1, 0) = u.y * v.x;
  result(1, 1) = u.y * v.y;
  result(1, 2) = u.y * v.z;

  result(2, 0) = u.z * v.x;
  result(2, 1) = u.z * v.y;
  result(2, 2) = u.z * v.z;

  return result;
}

float dotProd(const Vec3& u, const Vec3& v) {
  return u.x * v.x + u.y * v.y + u.z * v.z;
}

void Tree::insert(Node* node, const Particle& p) {
  lastInserted = node;
  if (node->children[0] != nullptr) {  // internal node
    insert(node->children[getChildId(node, p)], p);
    return;
  }
  // external node
  if (node->p == nullptr) {  // free spot
    node->p = &p;
    return;
  }
  // occupied spot
  createChildren(node);
  insert(node->children[getChildId(node, *node->p)], *node->p);
  insert(node->children[getChildId(node, p)], p);
  node->p = nullptr;
}

void createChildren(Node* parent) {
  for (int i = 0; i < 8; ++i) {
    float H = parent->H / 2;
    int x = (i >> 2) & 1;
    int y = (i >> 1) & 1;
    int z = i & 1;
    Vec3 low = parent->low + H * Vec3(float(x), float(y), float(z));
    parent->children[i] = new Node(low, H, parent);
  }
}

int getChildId(Node* root, const Particle& p) {
  Vec3 offset = 2 * (p.position - root->low) / root->H;
  int x = int(offset.x);
  int y = int(offset.y);
  int z = int(offset.z);
  int childId = (x << 2) | (y << 1) | z;
  return childId;
}

void DFSFree(Node* root) {
  if (root == nullptr) {
    return;
  }
  for (int i = 0; i < 8; ++i) {
    DFSFree(root->children[i]);
    delete root->children[i];
  }
}

void DFSPrint(Node* root) {
  if (root == nullptr) {
    return;
  }
  if (root->children[0] == nullptr) {
    std::cout << "[external node]\n";
    if (root->p) {
      Vec3 pos = root->p->position;
      std::cout << "particle position: " << pos.x << ' ' << pos.y << ' ' << pos.z << '\n';
    }
  } else {
    std::cout << "[internal node]\n";
    std::cout << "total node mass: " << root->M << '\n';
  }
  Vec3 low = root->low;
  std::cout << "low: " << low.x << ' ' << low.y << ' ' << low.z << '\n';
  std::cout << "H: " << root->H << '\n';
  std::cout << "ID: " << root->id << '\n';
  for (int i = 0; i < 8; ++i) {
    DFSPrint(root->children[i]);
  }
}

void calcQ(Node* node) {
  if (node->children[0] == nullptr) {
    node->Q = Mat3x3::zero();
    return;
  }
  auto Q = Mat3x3::zero();
  for (Node* c : node->children) {
    calcQ(c);
    Q += c->Q;
    Vec3 childCOM = c->p == nullptr ? c->COM : c->p->position;
    float childMass = c->p == nullptr ? c->M : c->p->mass;
    auto R = childCOM - node->COM;
    Q += (outerProd(R, R) * 3 - Mat3x3::I() * R.getMagnitudeSquared()) * childMass;
  }
  node->Q = Q;
}

void calcCOM(Node* node) {
  if (node == nullptr)
    return;

  // external node
  if (node->children[0] == nullptr) {
    if (node->p) {
      node->M = node->p->mass;
      node->COM = node->p->position;
    } else {
      node->M = 0;
      node->COM = Vec3::zero();
    }
    return;
  }

  // internal node
  node->M = 0;
  node->COM = Vec3::zero();
  for (int i = 0; i < 8; ++i) {
    calcCOM(node->children[i]);
    float m = node->children[i]->M;
    node->COM += m * node->children[i]->COM;
    node->M += m;
  }

  if (node->M > 0) {
    node->COM /= node->M;
  }
}

void findForce(const Node* root,
               Particle& p,
               float theta,
               float G,
               float eps,
               bool includeQuadrupole) {
  if (root->children[0] == nullptr) {  // external node
    if (root->p != nullptr && root->p != &p) {
      p.acceleration += gravity(root->p->position, root->p->mass, p.position, G, eps);
    }
    return;
  }
  // internal node
  if (root->H / dist(root->COM, p.position) < theta) {
    p.acceleration += gravity(root->COM, root->M, p.position, G, eps);
    if (includeQuadrupole) {
      auto rVec = p.position - root->COM;
      float r = std::sqrtf(rVec.getMagnitudeSquared());
      p.acceleration += G * (root->Q * rVec / std::powf(r, 5) -
                             (5 / 2.0f) * dotProd(rVec, root->Q * rVec) * rVec / std::powf(r, 7));
    }
    return;
  }
  for (const Node* c : root->children) {
    findForce(c, p, theta, G, eps, includeQuadrupole);
  }
}

float gravPotential(Vec3 r1, float m1, Vec3 r2, float m2, float G, float eps) {
  return -1 * G * m1 * m2 / std::sqrtf((r2 - r1).getMagnitudeSquared() + eps * eps);
}

void findPE(const Node* root,
            int i,
            const std::vector<Particle>& particles,
            float theta,
            float G,
            float eps,
            std::vector<float>& PEs,
            bool includeQuadrupole) {
  const Particle& p = particles[i];
  if (root->children[0] == nullptr) {  // external node
    if (root->p != nullptr && root->p != &p) {
      PEs[i] += gravPotential(root->p->position, root->p->mass, p.position, p.mass, G, eps);
    }
    return;
  }
  // internal node
  if (root->H / dist(root->COM, p.position) < theta) {
    PEs[i] += gravPotential(root->COM, root->M, p.position, p.mass, G, 0);
    if (includeQuadrupole) {
      auto rVec = p.position - root->COM;
      float r = rVec.getMagnitude();
      float quadPot = -0.5f * G / std::powf(r, 5) * dotProd(rVec, root->Q * rVec);
      PEs[i] += quadPot * p.mass;
    }
    return;
  }
  for (const Node* c : root->children) {
    findPE(c, i, particles, theta, G, eps, PEs, includeQuadrupole);
  }
}

Vec3 gravity(Vec3 r1, float m1, Vec3 r2, float G, float eps) {
  Vec3 r21 = r2 - r1;
  Vec3 g21 = -1 * G * m1 / std::powf(r21.getMagnitudeSquared() + eps * eps, 1.5f) * r21;
  return g21;
}

float dist(Vec3 a, Vec3 b) {
  return (a - b).getMagnitude();
}

void setHalfStepVelocities(std::vector<Particle>& particles, float dt) {
  std::for_each(std::execution::par_unseq, particles.begin(), particles.end(),
                [dt](Particle& p) { p.velocity += 0.5f * dt * p.acceleration; });
}

void setIntegerStepVelocities(std::vector<Particle>& particles, float dt) {
  std::for_each(std::execution::par, particles.begin(), particles.end(), [dt](Particle& p) {
    p.integerStepVelocity = p.velocity + 0.5f * dt * p.acceleration;
  });
}

void updateVelocities(std::vector<Particle>& particles, float dt) {
  std::for_each(std::execution::par_unseq, particles.begin(), particles.end(),
                [dt](Particle& p) { p.velocity += dt * p.acceleration; });
}

void updatePositions(std::vector<Particle>& particles, float dt) {
  std::for_each(std::execution::par_unseq, particles.begin(), particles.end(),
                [dt](Particle& p) { p.position += dt * p.velocity; });
}

BarnesHut::BarnesHut(const std::vector<Vec3>& state,
                     const std::vector<float>& masses,
                     std::function<Vec3(Vec3)> externalField,
                     std::function<float(Vec3)> externalPotential,
                     Vec3 low,
                     float H,
                     float G,
                     float softeningLength,
                     float theta,
                     bool includeQuadrupole,
                     bool useZOrdering)
    : N(int(masses.size())),
      externalField(externalField),
      externalPotential(externalPotential),
      PEContributions(N),
      low(low),
      H(H),
      G(G),
      eps(softeningLength),
      theta(theta),
      includeQuadrupole(includeQuadrupole),
      useZOrdering(useZOrdering) {
  for (int i = 0; i < N; ++i) {
    this->particles.emplace_back(state[i], state[N + i], masses[i]);
  }
}

BarnesHut::~BarnesHut() {}

void BarnesHut::run(StateRecorder& stateRecorder,
                    const int simLength,
                    float dt,
                    bool collectDiagnostics,
                    bool recordField) {
  SimInfo simInfo;
  simInfo.setInitialMomentum(particles);

  doStep();

  setHalfStepVelocities(particles, dt);

  for (int t = 0; t <= simLength; ++t) {
    std::cout << "progress: " << t << "/" << simLength << '\r';
    std::cout.flush();

    if (recordField) {
      stateRecorder.recordField(particles, 1, 1);
    }

    updatePositions(particles, dt);
    setIntegerStepVelocities(particles);

    stateRecorder.recordPositions(particles);
    auto expectedMomentum = simInfo.updateExpectedMomentum(totalExternalForce(), dt);
    stateRecorder.recordExpectedMomentum(expectedMomentum);
    stateRecorder.recordTotalMomentum(SimInfo::totalMomentum(particles));
    stateRecorder.recordTotalAngularMomentum(SimInfo::totalAngularMomentum(particles));
    float ke = SimInfo::kineticEnergy(particles);
    stateRecorder.recordEnergy(PE, ke);

    if (escapedBox()) {
      std::cout << "particle escaped box\n";
      break;
    }

    doStep();
    updateVelocities(particles, dt);
  }

  stateRecorder.flush();
}

double BarnesHut::getTreeConstructionSecs() {
  return treeConstructionSeconds;
}

void BarnesHut::clearAccelerations() {
  std::for_each(std::execution::par, particles.begin(), particles.end(),
                [](Particle& p) { p.acceleration = Vec3::zero(); });
}

void BarnesHut::calculateAccelerations(const Tree& tree) {
  std::for_each(std::execution::par, particles.begin(), particles.end(),
                [this, &tree](Particle& p) {
                  findForce(tree.getRoot(), p, theta, G, eps, includeQuadrupole);
                  p.acceleration += externalField(p.position);
                });
}

void BarnesHut::doStep() {
  clearAccelerations();
  auto start = std::chrono::high_resolution_clock::now();
  Tree tree(particles, low, H, useZOrdering);
  auto end = std::chrono::high_resolution_clock::now();
  treeConstructionSeconds += std::chrono::duration<double>(end - start).count();

  calcCOM(tree.root);
  if (includeQuadrupole) {
    calcQ(tree.root);
  }
  calculateAccelerations(tree);
  PE = calculatePE(tree);
}

bool isWithinBox(Vec3 pos, Vec3 low, float H) {
  return pos.x >= low.x && pos.x <= low.x + H && pos.y >= low.y && pos.y <= low.y + H &&
         pos.z >= low.z && pos.z <= low.z + H;
}

bool BarnesHut::escapedBox() {
  return std::any_of(particles.begin(), particles.end(),
                     [this](const Particle& p) { return !isWithinBox(p.position, low, H); });
}

Vec3 BH::BarnesHut::totalExternalForce() {
  return std::transform_reduce(
      std::execution::par, particles.begin(), particles.end(), Vec3::zero(), std::plus<>(),
      [this](const Particle& p) { return p.mass * externalField(p.position); });
}

float BarnesHut::calculatePE(const Tree& tree) {
  const auto nRange = std::views::iota(0, N);
  std::fill(std::execution::par, PEContributions.begin(), PEContributions.end(), 0.0f);
  std::for_each(std::execution::par, nRange.begin(), nRange.end(), [this, &tree](int i) {
    findPE(tree.getRoot(), i, particles, theta, G, eps, PEContributions, includeQuadrupole);
    Particle& p = particles[i];
    PEContributions[i] /= 2;
    PEContributions[i] += p.mass * externalPotential(p.position);
  });

  return std::reduce(std::execution::par, PEContributions.begin(), PEContributions.end());
}

Mat3x3 Mat3x3::zero() {
  return Mat3x3{std::array<float, 9>{0, 0, 0, 0, 0, 0, 0, 0, 0}};
}

Mat3x3 Mat3x3::I() {
  return Mat3x3{std::array<float, 9>{1, 0, 0, 0, 1, 0, 0, 0, 1}};
}

float& Mat3x3::operator()(int row, int col) {
  return data[row * 3 + col];
}

const float& Mat3x3::operator()(int row, int col) const {
  return data[row * 3 + col];
}

Mat3x3 Mat3x3::operator+(const Mat3x3& other) const {
  Mat3x3 result;
  for (int i = 0; i < 9; ++i)
    result.data[i] = data[i] + other.data[i];
  return result;
}

Mat3x3 Mat3x3::operator-(const Mat3x3& other) const {
  Mat3x3 result;
  for (int i = 0; i < 9; ++i)
    result.data[i] = data[i] - other.data[i];
  return result;
}

Mat3x3& Mat3x3::operator+=(const Mat3x3& other) {
  for (int i = 0; i < 9; ++i)
    data[i] += other.data[i];
  return *this;
}

Mat3x3 Mat3x3::operator*(float scalar) const {
  Mat3x3 result;
  for (int i = 0; i < 9; ++i)
    result.data[i] = data[i] * scalar;
  return result;
}

Mat3x3& Mat3x3::operator*=(float scalar) {
  for (int i = 0; i < 9; ++i)
    data[i] *= scalar;
  return *this;
}

}  // namespace BH

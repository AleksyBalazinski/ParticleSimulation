#include "barnesHut.h"
#include <algorithm>
#include <cmath>
#include <execution>
#include <iostream>
#include "simInfo.h"
#include "stateRecorder.h"

namespace BH {

int nodeID = 0;

Node::Node(Vec3 low, float H, Particle* p)
    : low(low), H(H), p(p), COM(), M(0.0f), children({nullptr}) {
  id = nodeID++;
}

Node::~Node() {}

Tree::Tree(const std::vector<Particle>& particles, Vec3 low, float H) {
  root = new Node(low, H);
  createChildren(root->children, low, H);

  for (const auto& p : particles) {
    insert(root, p);
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

void insert(Node* root, const Particle& p) {
  if (root->children[0] != nullptr) {  // internal node
    updateInternalNode(root, p);
    insert(root->children[getChildId(root, p)], p);
    return;
  }
  // external node
  if (root->p == nullptr) {  // free spot
    root->p = &p;
    return;
  }
  // occupied spot
  createChildren(root->children, root->low, root->H);
  updateInternalNode(root, *root->p);
  updateInternalNode(root, p);
  insert(root->children[getChildId(root, *root->p)], *root->p);
  insert(root->children[getChildId(root, p)], p);
  root->p = nullptr;
}

void createChildren(std::array<Node*, 8>& children, Vec3 parentLow, float parentH) {
  for (int i = 0; i < 8; ++i) {
    float H = parentH / 2;
    int x = (i >> 2) & 1;
    int y = (i >> 1) & 1;
    int z = i & 1;
    Vec3 low = parentLow + H * Vec3(float(x), float(y), float(z));
    children[i] = new Node(low, H);
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

void updateInternalNode(Node* node, const Particle& p) {
  float totalMass = node->M + p.mass;
  node->COM = (node->M * node->COM + p.mass * p.position) / totalMass;
  node->M = totalMass;
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

void findForce(const Node* root, Particle& p, float theta, float G, float eps) {
  if (root->children[0] == nullptr) {  // external node
    if (root->p != nullptr && root->p != &p) {
      p.acceleration += gravity(root->p->position, root->p->mass, p.position, G, eps);
    }
    return;
  }
  // internal node
  if (root->H / dist(root->COM, p.position) < theta) {
    p.acceleration += gravity(root->COM, root->M, p.position, G, eps);
    return;
  }
  for (const Node* c : root->children) {
    findForce(c, p, theta, G, eps);
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
                     float theta)
    : N(int(masses.size())),
      externalField(externalField),
      externalPotential(externalPotential),
      low(low),
      H(H),
      G(G),
      eps(softeningLength),
      theta(theta) {
  for (int i = 0; i < N; ++i) {
    this->particles.emplace_back(state[i], state[N + i], masses[i]);
  }
}

BarnesHut::~BarnesHut() {}

void BarnesHut::run(int simLength, float dt) {
  StateRecorder stateRecorder("output.dat", N, simLength + 1, "energy.txt", "momentum.txt",
                              "expected_momentum.txt", "angular_momentum.txt", "");
  SimInfo simInfo;
  simInfo.setInitialMomentum(particles);

  doStep();
  setHalfStepVelocities(particles, dt);

  for (int t = 0; t <= simLength; ++t) {
    std::cout << "progress: " << float(t) / simLength << '\r';
    std::cout.flush();

    updatePositions(particles, dt);
    setIntegerStepVelocities(particles);

    stateRecorder.recordPositions(particles);
    auto expectedMomentum = simInfo.updateExpectedMomentum(totalExternalForce(), dt);
    stateRecorder.recordExpectedMomentum(expectedMomentum);
    stateRecorder.recordTotalMomentum(SimInfo::totalMomentum(particles));
    stateRecorder.recordTotalAngularMomentum(SimInfo::totalAngularMomentum(particles));
    float pe = 0.0f;  // TODO n^2 complexity for naive computation -> borrow a grid from PM or
    // approximate from tree ??
    float ke = SimInfo::kineticEnergy(particles);
    stateRecorder.recordEnergy(pe, ke);

    if (escapedBox()) {
      std::cout << "particle escaped box\n";
      break;
    }

    doStep();
    updateVelocities(particles, dt);
  }

  stateRecorder.flush();
}

void BarnesHut::clearAccelerations() {
  std::for_each(std::execution::par, particles.begin(), particles.end(),
                [](Particle& p) { p.acceleration = Vec3::zero(); });
}

void BarnesHut::calculateAccelerations(const Tree& tree) {
  std::for_each(std::execution::par, particles.begin(), particles.end(),
                [this, &tree](Particle& p) {
                  findForce(tree.getRoot(), p, theta, G, eps);
                  p.acceleration += externalField(p.position);
                });
}

void BarnesHut::doStep() {
  clearAccelerations();
  Tree tree(particles, low, H);
  calculateAccelerations(tree);
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

}  // namespace BH

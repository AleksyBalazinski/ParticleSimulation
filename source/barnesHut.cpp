#include "barnesHut.h"
#include <algorithm>
#include <cmath>
#include <execution>
#include <iostream>
#include "stateRecorder.h"

namespace BH {

int nodeID = 0;

// Particle::Particle(Vec3 position, Vec3 velocity, float mass)
//     : position(position), velocity(velocity), mass(mass) {}

Node::Node(Vec3 low, float H, Particle* p) : low(low), H(H), p(p), COM(), M(0.0f) {
  id = nodeID++;
  for (int i = 0; i < 8; ++i) {
    children[i] = nullptr;
  }
}

Node::~Node() {
  // std::cout << "removing " << id << '\n';
}

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
  node->COM = node->M * node->COM + p.mass * p.position;
  node->M += p.mass;
  node->COM /= node->M;
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

void findForce(const Node* root, Particle& p, float threshold, float G) {
  if (root->children[0] == nullptr) {  // external node
    if (root->p != nullptr && root->p != &p) {
      p.acceleration += gravity(root->p->position, root->p->mass, p.position, p.mass, G) / p.mass;
    }
    return;
  }
  // internal node
  if (dist(root->COM, p.position) / root->H > threshold) {
    p.acceleration += gravity(root->COM, root->M, p.position, p.mass, G) / p.mass;
    return;
  }
  for (const Node* c : root->children) {
    findForce(c, p, threshold, G);
  }
}

Vec3 gravity(Vec3 r1, float m1, Vec3 r2, float m2, float G) {
  Vec3 r21 = r2 - r1;
  return -1 * G * m1 * m2 / std::powf(r21.getMagnitude() + 0.01f, 3.0f) * r21;
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
                     Vec3 low,
                     float H,
                     float G)
    : N(int(masses.size())), low(low), H(H), G(G) {
  for (int i = 0; i < N; ++i) {
    this->particles.emplace_back(state[i], state[N + i], masses[i]);
  }
}

BarnesHut::~BarnesHut() {}

void BarnesHut::run(int simLength, float dt) {
  bool withinBox = true;
  StateRecorder stateRecorder("output.dat", N, simLength + 1, "", "", "", "", "");

  setHalfStepVelocities(particles, dt);
  doStep();
  for (int t = 0; t <= simLength; ++t) {
    std::cout << "progress: " << float(t) / simLength << '\r';
    std::cout.flush();

    updatePositions(particles, dt);
    stateRecorder.recordPositions(particles);
    if (escapedBox()) {
      withinBox = false;
      break;
    }

    doStep();
    updateVelocities(particles, dt);
  }
  if (!withinBox) {
    std::cout << "particle escaped box\n";
  } else {
    stateRecorder.flush();
  }
}

void BarnesHut::clearAccelerations() {
  for (auto& p : particles) {
    p.acceleration = Vec3();
  }
}

void BarnesHut::calculateAccelerations(const Tree& tree) {
  for (auto& p : particles) {
    findForce(tree.getRoot(), p, 1, G);  // TODO variable threshold
  }
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

}  // namespace BH

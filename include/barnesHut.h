#pragma once

#include <array>
#include <vector>
#include "particle.h"
#include "vec3.h"

namespace BH {

// struct Particle {
//   Vec3 position;
//   Vec3 velocity;
//   Vec3 acceleration;
//   float mass;
//   Vec3 integerStepVelocity;

//   Particle(Vec3 position, Vec3 velocity, float mass);
// };

struct Node {
  Node(Vec3 low, float H, Particle* p = nullptr);
  ~Node();
  Vec3 COM;
  float M;
  Vec3 low;
  float H;
  std::array<Node*, 8> children;
  const Particle* p;
  int id;  // TODO remove
};

void insert(Node* root, const Particle& p);
void createChildren(std::array<Node*, 8>& children, Vec3 parentLow, float parentH);
int getChildId(Node* root, const Particle& p);
void updateInternalNode(Node* node, const Particle& p);
void DFSFree(Node* root);
void DFSPrint(Node* root);

void findForce(const Node* root, Particle& p, float threshold, float G);
Vec3 gravity(Vec3 r1, float m1, Vec3 r2, float m2, float G);
float dist(Vec3 a, Vec3 b);

void setHalfStepVelocities(std::vector<Particle>& particles, float dt = 1.0f);
void setIntegerStepVelocities(std::vector<Particle>& particles, float dt = 1.0f);
void updateVelocities(std::vector<Particle>& particles, float dt = 1.0f);
void updatePositions(std::vector<Particle>& particles, float dt = 1.0f);

class Tree {
 public:
  Tree(const std::vector<Particle>& particles, Vec3 low, float H);
  ~Tree();
  void print();
  const Node* getRoot() const;

 private:
  Node* root;
};

class BarnesHut {
 public:
  BarnesHut(const std::vector<Vec3>& state,
            const std::vector<float>& masses,
            Vec3 low,
            float H,
            float G);
  ~BarnesHut();

  void run(int simLength, float dt);

 private:
  void clearAccelerations();
  void calculateAccelerations(const Tree& t);
  void doStep();

  bool escapedBox();

  int N;
  std::vector<Particle> particles;
  Vec3 low;
  float H;
  float G;
};

}  // namespace BH

#pragma once

#include <array>
#include <functional>
#include <vector>
#include "particle.h"
#include "stateRecorder.h"
#include "vec3.h"

namespace BH {

struct Mat3x3 {
  static Mat3x3 zero();

  static Mat3x3 I();

  float& operator()(int row, int col);

  const float& operator()(int row, int col) const;

  Mat3x3 operator+(const Mat3x3& other) const;

  Mat3x3 operator-(const Mat3x3& other) const;

  Mat3x3& operator+=(const Mat3x3& other);

  Mat3x3 operator*(float scalar) const;

  Mat3x3& operator*=(float scalar);

  std::array<float, 9> data;
};

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
  Mat3x3 Q;
};

Vec3 operator*(const Mat3x3& m, const Vec3& v);
Mat3x3 outerProd(const Vec3& u, const Vec3& v);
float dotProd(const Vec3& u, const Vec3& v);

void insert(Node* root, const Particle& p);
void createChildren(std::array<Node*, 8>& children, Vec3 parentLow, float parentH);
int getChildId(Node* root, const Particle& p);
void updateInternalNode(Node* node, const Particle& p);
void DFSFree(Node* root);
void DFSPrint(Node* root);

void calcQ(Node* node);
void findForce(const Node* root,
               Particle& p,
               float theta,
               float G,
               float eps,
               bool includeQuadrupole);
Vec3 gravity(Vec3 r1, float m1, Vec3 r2, float G, float eps);
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

  Node* root;
};

class BarnesHut {
 public:
  BarnesHut(const std::vector<Vec3>& state,
            const std::vector<float>& masses,
            std::function<Vec3(Vec3)> externalField,
            std::function<float(Vec3)> externalPotential,
            Vec3 low,
            float H,
            float G,
            float softeningLength,
            float theta,
            bool includeQuadrupole = false);
  ~BarnesHut();

  void run(StateRecorder& stateRecorder,
           const int simLength,
           float dt,
           bool collectDiagnostics = false,
           bool recordField = false);

 private:
  void clearAccelerations();
  void calculateAccelerations(const Tree& t);
  void doStep();

  bool escapedBox();

  Vec3 totalExternalForce();
  float calculatePE(const Tree& t);

  int N;
  std::vector<Particle> particles;

  std::function<Vec3(Vec3)> externalField;
  std::function<float(Vec3)> externalPotential;
  std::vector<float> PEContributions;
  float PE;

  Vec3 low;
  float H;
  float G;
  float eps;
  float theta;
  bool includeQuadrupole;
};

}  // namespace BH

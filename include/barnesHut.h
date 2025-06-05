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
  Node(Vec3 low, float H, Node* parent, Particle* p = nullptr);
  ~Node();
  Vec3 COM;
  float M;
  Vec3 low;
  float H;
  std::array<Node*, 8> children;
  Node* parent;
  const Particle* p;
  int id;  // TODO remove
  Mat3x3 Q;
};

Vec3 operator*(const Mat3x3& m, const Vec3& v);
Mat3x3 outerProd(const Vec3& u, const Vec3& v);
float dotProd(const Vec3& u, const Vec3& v);

void createChildren(Node* parent);
int getChildId(Node* root, const Particle& p);
void DFSFree(Node* root);
void DFSPrint(Node* root);

void calcQ(Node* node);
void calcCOM(Node* node);
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
  Tree(const std::vector<Particle>& particles, Vec3 low, float H, bool useZOrdering);
  ~Tree();
  void print();
  const Node* getRoot() const;
  void initZCodes(const std::vector<Particle>& particles);
  void insert(Node* root, const Particle& p);

  Vec3 compBoxLow;
  float compBoxSize;
  std::vector<uint32_t> zCodes;
  Node* root;
  Node* lastInserted;
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
            bool includeQuadrupole = false,
            bool useZOrdering = false);
  ~BarnesHut();

  void run(StateRecorder& stateRecorder,
           const int simLength,
           float dt,
           bool collectDiagnostics = false,
           bool recordField = false);

  double getTreeConstructionSecs();

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
  bool useZOrdering;
  double treeConstructionSeconds{};
};

}  // namespace BH

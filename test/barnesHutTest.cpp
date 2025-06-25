#include <gtest/gtest.h>
#include "barnesHut.h"

using namespace BH;

TEST(barnesHut, tree) {
  std::vector<Particle> particles;
  Particle p1(Vec3(1, 1, 1), Vec3(), 1);
  Particle p2(Vec3(1, 1, 3), Vec3(), 2);
  Particle p3(Vec3(9, 9, 9), Vec3(), 3);
  particles.push_back(p1);
  particles.push_back(p2);
  particles.push_back(p3);

  Vec3 low(0, 0, 0);
  float H = 10;
  Tree tree(particles, low, H, false);
  tree.print();
}
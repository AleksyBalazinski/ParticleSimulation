#pragma once

#include <vector>
#include "vec3.h"

class Grid {
 private:
  int gridPoints;
  int length;
  std::vector<double> density;
  std::vector<Vec3> field;

 public:
  Grid(int gridPoints);

  void assignDensity(int x, int y, int z, double density);

  void clearDensity();

  void assignField(int x, int y, int z, Vec3 fieldVal);

  Vec3 getField(int x, int y, int z);

  int getLength() const { return length; }

  int getGridPoints() const { return gridPoints; }

  const std::vector<double>& getDensity() const { return density; }
};

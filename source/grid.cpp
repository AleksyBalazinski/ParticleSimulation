#include "grid.h"
#include <algorithm>

Grid::Grid(int gridPoints)
    : gridPoints(gridPoints),
      length(gridPoints * gridPoints * gridPoints),
      density(length),
      field(length) {
  for (int i = 0; i < density.size(); i++) {
    density[i] = 0.0;
  }
}

void Grid::assignDensity(int x, int y, int z, double d) {
  density[x * gridPoints * gridPoints + y * gridPoints + z] += d;
}

void Grid::clearDensity() {
  std::fill(density.begin(), density.end(), 0);
}

void Grid::assignField(int x, int y, int z, Vec3 fieldVal) {
  field[x * gridPoints * gridPoints + y * gridPoints + z] = fieldVal;
}

Vec3 Grid::getField(int x, int y, int z) {
  return field[x * gridPoints * gridPoints + y * gridPoints + z];
}

#pragma once

#include <functional>
#include <vector>
#include "grid.h"
#include "stateRecorder.h"
#include "vec3.h"

enum class InterpolationScheme { NGP, CIC };
enum class FiniteDiffScheme { TWO_POINT, FOUR_POINT };

class PMMethod {
 private:
  Grid& grid;
  std::vector<Vec3>& state;
  std::vector<double>& masses;
  std::function<Vec3(Vec3)> externalField;
  double H;
  double DT;
  double G;
  InterpolationScheme is;
  FiniteDiffScheme fds;

  void reassignDensity();
  Vec3 getField(double x, double y, double z);

  void updateAccelerations(std::vector<Vec3>& accelerations, StateRecorder& sr);

 public:
  PMMethod(std::vector<Vec3>& state,
           std::vector<double>& masses,
           std::function<Vec3(Vec3)> externalField,
           double H,
           double DT,
           double G,
           InterpolationScheme is,
           FiniteDiffScheme fds,
           Grid& grid);

  std::string run(const int simLength,
                  bool collectDiagnostics = false,
                  const char* positionsPath = "output.txt",
                  const char* energyPath = "energy.txt",
                  const char* momentumPath = "momentum.txt");
};

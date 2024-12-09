#include "unit_conversions.h"
#include <algorithm>
#include <numbers>

Vec3 positionToCodeUntits(const Vec3& pos, double H) {
  return pos / H;
}

Vec3 velocityToCodeUntits(const Vec3& v, double H, double DT) {
  return DT * v / H;
}

double densityToCodeUnits(double density, double DT, double G) {
  return DT * DT * 4 * std::numbers::pi * G * density;
}

void stateToCodeUnits(std::vector<Vec3>& state, double H, double DT) {
  int N = (int)state.size() / 2;
  std::transform(state.begin(), state.begin() + N, state.begin(),
                 [H](const Vec3& pos) { return positionToCodeUntits(pos, H); });
  std::transform(state.begin() + N, state.end(), state.begin() + N,
                 [H, DT](const Vec3& v) { return velocityToCodeUntits(v, H, DT); });
}

void velocitiesToCodeUnits(std::vector<Vec3>& velocities, double H, double DT) {
  std::transform(velocities.begin(), velocities.end(), velocities.begin(),
                 [H, DT](const Vec3& v) { return velocityToCodeUntits(v, H, DT); });
}

Vec3 positionToOriginalUnits(const Vec3& pos, double H) {
  return H * pos;
}

Vec3 velocityToOriginalUnits(const Vec3& v, double H, double DT) {
  return H * v / DT;
}

void stateToOriginalUnits(std::vector<Vec3>& state, double H, double DT) {
  int N = (int)state.size() / 2;
  std::transform(state.begin(), state.begin() + N, state.begin(),
                 [H](const Vec3& pos) { return positionToOriginalUnits(pos, H); });
  std::transform(state.begin() + N, state.end(), state.begin() + N,
                 [H, DT](const Vec3& v) { return velocityToOriginalUnits(v, H, DT); });
}

void velocitiesToOriginalUnits(std::vector<Vec3>& velocities, double H, double DT) {
  std::transform(velocities.begin(), velocities.end(), velocities.begin(),
                 [H, DT](const Vec3& v) { return velocityToOriginalUnits(v, H, DT); });
}
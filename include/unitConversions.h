#pragma once

#include <numbers>
#include <vector>
#include "particle.h"
#include "vec3.h"

inline Vec3 positionToCodeUntits(const Vec3& pos, float H) {
  return pos / H;
}
inline Vec3 positionToOriginalUnits(const Vec3& pos, float H) {
  return H * pos;
}

inline Vec3 velocityToCodeUntits(const Vec3& v, float H, float DT) {
  return DT * v / H;
}
inline Vec3 velocityToOriginalUnits(const Vec3& v, float H, float DT) {
  return H * v / DT;
}

inline Vec3 accelerationToCodeUnits(const Vec3& a, float H, float DT) {
  return DT * DT * a / H;
}

inline Vec3 accelerationToOriginalUnits(const Vec3& a, float H, float DT) {
  return H * a / (DT * DT);
}

inline float densityToCodeUnits(float density, float DT, float G) {
  return DT * DT * 4 * std::numbers::pi_v<float> * G * density;
}
inline float densityToOriginalUnits(float density, float DT, float G) {
  return density / (DT * DT * 4 * std::numbers::pi_v<float> * G);
}

inline float potentialToOriginalUnits(float potential, float H, float DT) {
  return potential * H * H / (DT * DT);
}

inline float massToCodeUnits(float m, float H, float DT, float G) {
  return DT * DT * 4 * std::numbers::pi_v<float> * G / (H * H * H) * m;
}
inline float massToOriginalUnits(float m, float H, float DT, float G) {
  return (H * H * H) / (DT * DT * 4 * std::numbers::pi_v<float> * G) * m;
}

inline float lengthToCodeUnits(float x, float H) {
  return x / H;
}

void stateToCodeUnits(std::vector<Vec3>& state, float H, float DT);
void stateToOriginalUnits(std::vector<Vec3>& state, float H, float DT);

void stateToCodeUnits(std::vector<Particle>& particles, float H, float DT);
void stateToOriginalUnits(std::vector<Particle>& particles, float H, float DT);

void velocitiesToCodeUnits(std::vector<Vec3>& velocities, float H, float DT);
void velocitiesToOriginalUnits(std::vector<Vec3>& velocities, float H, float DT);

void integerStepVelocitiesToOriginalUnits(std::vector<Particle>& particles, float H, float DT);
void integerStepVelocitiesToCodeUnits(std::vector<Particle>& particles, float H, float DT);

void massToCodeUnits(std::vector<Particle>& particles, float H, float DT, float G);
void massToOriginalUnits(std::vector<Particle>& particles, float H, float DT, float G);

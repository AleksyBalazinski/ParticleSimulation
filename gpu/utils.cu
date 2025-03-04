#include <cmath>
#include "utils.cuh"

__device__ Vec3 externalFieldBulge(Vec3 pos, Vec3 bulge, float rb, float mb, float G) {
  Vec3 d(pos.x - bulge.x, pos.y - bulge.y, pos.z - bulge.z);
  float r = norm3df(d.x, d.y, d.z);
  Vec3 dir(d.x / r, d.y / r, d.z / r);
  float g;
  if (r > rb) {
    g = -G * mb / (r * r);
  } else {
    g = -(G * mb / (rb * rb * rb)) * r * (4 - 3 * r / rb);
  }

  return Vec3(g * dir.x, g * dir.y, g * dir.z);
}

__host__ Vec3 externalFieldBulgeHost(Vec3 pos, Vec3 bulge, float rb, float mb, float G) {
  float r = (pos - bulge).getMagnitude();
  Vec3 dir = (pos - bulge) / r;
  float g;
  if (r > rb) {
    g = -G * mb / (r * r);
  } else {
    g = -(G * mb / std::powf(rb, 3)) * r * (4 - 3 * r / rb);
  }

  return g * dir;
}

__host__ float newton(std::function<float(float)> f,
                      std::function<float(float)> df,
                      float guess,
                      int maxIter,
                      float tolerance) {
  int i = 0;
  float err = tolerance + 1;
  float x = guess;
  float xPrev;
  while (err > tolerance && i < maxIter) {
    ++i;
    xPrev = x;
    x = x - f(x) / df(x);
    err = std::abs(x - xPrev);
  }

  return x;
}

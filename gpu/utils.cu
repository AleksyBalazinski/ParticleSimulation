#include "utils.cuh"

__host__ __device__ Vec3 externalFieldBulge(Vec3 pos, Vec3 bulge, float rb, float mb, float G) {
  float r = (pos - bulge).getMagnitude();
  Vec3 dir = (pos - bulge) / r;
  float g;
  if (r > rb) {
    g = -G * mb / (r * r);
  } else {
    g = -(G * mb / powf(rb, 3)) * r * (4 - 3 * r / rb);
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

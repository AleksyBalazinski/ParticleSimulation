#include <cmath>
#include "utils.cuh"

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

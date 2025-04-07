#pragma once

#include <cufft.h>

__device__ inline cufftComplex mul(cufftComplex a, cufftComplex b) {
  return make_cuComplex(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}
#pragma once

#include <cufft.h>
#include <tuple>
#include "FFTAdapter.h"

class CuFFTAdapter {
 public:
  CuFFTAdapter(std::tuple<int, int, int> dims);
  void free();

  void fft(cuComplex* d_in, cuComplex* d_out);
  void ifft(cuComplex* d_in, cuComplex* d_out);

 private:
  cufftHandle plan;
  int length;
};
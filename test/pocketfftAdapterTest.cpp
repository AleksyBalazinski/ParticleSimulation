#include <gtest/gtest.h>
#include <complex>
#include <vector>
#include "kissFFTAdapter.h"
#include "pocketfftAdapter.h"

TEST(pocketfftAdapter, forwardAndBackward) {
  int dims[] = {2, 2, 2};
  int ndim = 3;
  std::vector<std::complex<float>> original = {1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<std::complex<float>> transformed(8, 0.0f);
  PocketfftAdapter<float> adapter(dims, ndim);
  adapter.fft(original, transformed);
  for (auto& e : transformed) {
    std::cerr << e << ' ';
  }
  std::cerr << '\n';
  std::vector<std::complex<float>> invTransformed(8, 0.0f);
  adapter.ifft(transformed, invTransformed);
  for (auto& e : invTransformed) {
    std::cerr << e << ' ';
  }
}

TEST(kissFFTApapter, forwardAndBackward) {
  int dims[] = {2, 2, 2};
  int ndim = 3;
  std::vector<std::complex<float>> original = {1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<std::complex<float>> transformed(8, 0.0f);
  KissFFTAdapter<float> adapter(dims, ndim);
  adapter.fft(original, transformed);
  for (auto& e : transformed) {
    std::cerr << e << ' ';
  }
  std::cerr << '\n';
  std::vector<std::complex<float>> invTransformed(8, 0.0f);
  adapter.ifft(transformed, invTransformed);
  for (auto& e : invTransformed) {
    std::cerr << e << ' ';
  }
}

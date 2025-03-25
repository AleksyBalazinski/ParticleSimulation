#include <gtest/gtest.h>
#include <complex>
#include <vector>
#include "kissFFTAdapter.h"

TEST(kissFFTApapter, forwardAndBackward) {
  int dims[] = {2, 2, 2};
  int ndim = 3;
  std::vector<std::complex<float>> original = {1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<std::complex<float>> transformed(8, 0.0f);
  KissFFTAdapter<float> adapter(dims, ndim);

  adapter.fft(original, transformed);

  std::vector<std::complex<float>> invTransformed(8, 0.0f);
  adapter.ifft(transformed, invTransformed);

  for (int i = 0; i < original.size(); ++i) {
    ASSERT_NEAR(invTransformed[i].real(), original[i].real(), 1e-6f);
    ASSERT_NEAR(invTransformed[i].imag(), original[i].imag(), 1e-6f);
  }
}

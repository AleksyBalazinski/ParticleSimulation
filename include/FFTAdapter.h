#pragma once

#include <complex>
#include <vector>

template <typename T>
class FFTAdapter {
  using Ty = std::complex<T>;

 public:
  virtual ~FFTAdapter() {};

  virtual std::vector<std::complex<T>>& fft(std::vector<std::complex<T>>& in,
                                            std::vector<std::complex<T>>& out) = 0;
  virtual std::vector<std::complex<T>>& ifft(std::vector<std::complex<T>>& in,
                                             std::vector<std::complex<T>>& out) = 0;
};
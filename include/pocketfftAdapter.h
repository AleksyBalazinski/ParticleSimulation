#pragma once

#include <pocketfft_hdronly.h>
#include <algorithm>
#include <execution>
#include "FFTAdapter.h"

template <typename T>
class PocketfftAdapter : public FFTAdapter<T> {
  using CpxT = std::complex<T>;

 private:
  pocketfft::shape_t shape;
  pocketfft::stride_t stride;
  pocketfft::shape_t axes;
  int length;

 public:
  PocketfftAdapter(int* dims, int ndim)  // TODO generalize to other than 3D
      : shape(dims, dims + ndim),
        stride{sizeof(CpxT), dims[0] * sizeof(CpxT), dims[0] * dims[1] * sizeof(CpxT)},
        axes{0, 1, 2} {
    length = 1;
    for (int i = 0; i < ndim; i++) {
      length *= dims[i];
    }
  }

  virtual std::vector<CpxT>& fft(std::vector<CpxT>& in, std::vector<CpxT>& out) override {
    pocketfft::c2c(shape, stride, stride, axes, /*forward*/ true, in.data(), out.data(),
                   static_cast<T>(1.0));
    return out;
  }

  virtual std::vector<CpxT>& ifft(std::vector<CpxT>& in, std::vector<CpxT>& out) override {
    pocketfft::c2c(shape, stride, stride, axes, /*forward*/ false, in.data(), out.data(),
                   static_cast<T>(1.0));

    std::for_each(std::execution::par_unseq, out.begin(), out.end(),
                  [length = static_cast<T>(length)](CpxT& x) { x /= length; });
    return out;
  }
};

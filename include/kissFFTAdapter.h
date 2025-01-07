#pragma once

#include <kiss_fftnd.h>
#include "FFTAdapter.h"

template <typename T>
class KissFFTAdapter : public FFTAdapter<T> {
  using CpxT = std::complex<T>;

 private:
  kiss_fftnd_cfg cfg;
  kiss_fftnd_cfg cfgInv;

 public:
  KissFFTAdapter(int* dims, int ndims);
  ~KissFFTAdapter();

  virtual std::vector<CpxT>& fft(std::vector<CpxT>& in, std::vector<CpxT>& out) override {
    kiss_fftnd(cfg, reinterpret_cast<kiss_fft_cpx*>(in.data()),
               reinterpret_cast<kiss_fft_cpx*>(out.data()));
    return out;
  }

  virtual std::vector<CpxT>& ifft(std::vector<CpxT>& in, std::vector<CpxT>& out) override {
    kiss_fftnd(cfgInv, reinterpret_cast<kiss_fft_cpx*>(in.data()),
               reinterpret_cast<kiss_fft_cpx*>(out.data()));
    return out;
  }
};

template <typename T>
KissFFTAdapter<T>::KissFFTAdapter(int* dims, int ndims)
    : cfg(kiss_fftnd_alloc(dims, ndims, false, nullptr, nullptr)),
      cfgInv(kiss_fftnd_alloc(dims, ndims, true, nullptr, nullptr)) {}

template <typename T>
KissFFTAdapter<T>::~KissFFTAdapter() {
  kiss_fft_free(cfg);
  kiss_fft_free(cfgInv);
}

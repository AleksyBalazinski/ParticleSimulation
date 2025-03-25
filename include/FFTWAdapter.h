#pragma once

#include <fftw3.h>
#include <array>
#include "FFTAdapter.h"

class FFTWAdapter : public FFTAdapter<float> {
  using CpxT = std::complex<float>;

 public:
  FFTWAdapter(std::array<int, 3> dims);

  ~FFTWAdapter();

  virtual std::vector<CpxT>& fft(std::vector<CpxT>& in, std::vector<CpxT>& out) override;

  virtual std::vector<CpxT>& ifft(std::vector<CpxT>& in, std::vector<CpxT>& out) override;

 private:
  void copyToInputBuf(const std::vector<CpxT>& in);

  void copyFromOutputBuf(std::vector<CpxT>& out);

  int length;
  fftwf_complex* inBuf;
  fftwf_complex* outBuf;
  fftwf_plan fwdPlan;
  fftwf_plan invPlan;
};
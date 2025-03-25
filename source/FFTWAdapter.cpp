#include "FFTWAdapter.h"

#include <algorithm>
#include <execution>
#include <numeric>
#include <ranges>

using CpxT = std::complex<float>;

FFTWAdapter::FFTWAdapter(std::array<int, 3> dims) {
  length = std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<int>());
  inBuf = (fftwf_complex*)fftwf_malloc(length * sizeof(fftwf_complex));
  outBuf = (fftwf_complex*)fftwf_malloc(length * sizeof(fftwf_complex));

  fwdPlan = fftwf_plan_dft_3d(dims[0], dims[1], dims[2], inBuf, outBuf, FFTW_FORWARD, FFTW_MEASURE);
  invPlan =
      fftwf_plan_dft_3d(dims[0], dims[1], dims[2], inBuf, outBuf, FFTW_BACKWARD, FFTW_MEASURE);
}

FFTWAdapter::~FFTWAdapter() {
  fftwf_destroy_plan(fwdPlan);
  fftwf_destroy_plan(invPlan);
  fftwf_free(inBuf);
  fftwf_free(outBuf);
}

std::vector<CpxT>& FFTWAdapter::fft(std::vector<CpxT>& in, std::vector<CpxT>& out) {
  copyToInputBuf(in);
  fftwf_execute(fwdPlan);
  copyFromOutputBuf(out);
  return out;
}

std::vector<CpxT>& FFTWAdapter::ifft(std::vector<CpxT>& in, std::vector<CpxT>& out) {
  copyToInputBuf(in);
  fftwf_execute(invPlan);
  copyFromOutputBuf(out);

  std::for_each(std::execution::par_unseq, out.begin(), out.end(),
                [this](CpxT& x) { x /= float(length); });

  return out;
}

void FFTWAdapter::copyToInputBuf(const std::vector<CpxT>& in) {
  auto inRange = std::ranges::views::iota(0, int(in.size()));
  std::for_each(std::execution::par_unseq, inRange.begin(), inRange.end(), [this, &in](int i) {
    inBuf[i][0] = in[i].real();
    inBuf[i][1] = in[i].imag();
  });
}

void FFTWAdapter::copyFromOutputBuf(std::vector<CpxT>& out) {
  auto outRange = std::ranges::views::iota(0, int(out.size()));
  std::for_each(std::execution::par_unseq, outRange.begin(), outRange.end(),
                [this, &out](int i) { out[i] = {outBuf[i][0], outBuf[i][1]}; });
}

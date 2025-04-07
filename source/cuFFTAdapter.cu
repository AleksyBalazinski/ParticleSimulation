#include <numeric>
#include "cuFFTAdapter.cuh"

#define numBlocks(n) (((BLOCK_SIZE) + (n) - 1) / (BLOCK_SIZE))
#define BLOCK_SIZE 1024

CuFFTAdapter::CuFFTAdapter(std::tuple<int, int, int> dims)
    : length(std::get<0>(dims) * std::get<1>(dims) * std::get<2>(dims)) {
  cufftPlan3d(&plan, std::get<2>(dims), std::get<1>(dims), std::get<0>(dims), CUFFT_C2C);
}

void CuFFTAdapter::free() {
  cufftDestroy(plan);
}

void CuFFTAdapter::fft(cuComplex* d_in, cuComplex* d_out) {
  cufftExecC2C(plan, d_in, d_out, CUFFT_FORWARD);
}

__global__ void scaleAfterInverse(cuComplex* data, int length) {
  for (int idx = threadIdx.x + blockDim.x * blockIdx.x; idx < length;
       idx += gridDim.x * blockDim.x) {
    data[idx].x /= length;
    data[idx].y /= length;
  }
}

void CuFFTAdapter::ifft(cuComplex* d_in, cuComplex* d_out) {
  cufftExecC2C(plan, d_in, d_out, CUFFT_INVERSE);
  scaleAfterInverse<<<numBlocks(length), BLOCK_SIZE>>>(d_out, length);
}

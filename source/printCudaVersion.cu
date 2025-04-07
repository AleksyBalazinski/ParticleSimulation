#include "gpu.h"

#include <iostream>

void printCudaVersion() {
  std::cout << "---------- CUDA INFO ----------\n";
  std::cout << "CUDA Compiled version: " << __CUDACC_VER_MAJOR__ << '.' << __CUDACC_VER_MINOR__
            << '.' << __CUDACC_VER_BUILD__ << '\n';

  int runtime_ver;
  cudaRuntimeGetVersion(&runtime_ver);
  std::cout << "CUDA Runtime version: " << runtime_ver << '\n';

  int driver_ver;
  cudaDriverGetVersion(&driver_ver);
  std::cout << "CUDA Driver version: " << driver_ver << '\n';
  std::cout << "-------------------------------\n\n";
}
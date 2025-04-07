#include <iostream>
#include "demos.h"
#ifdef CUDA
#include "gpu.h"
#endif
int main() {
#ifdef CUDA
  printCudaVersion();
#endif
  try {
    // smallSimPP();
    // bigSimPP();
    // smallSimP3M();
    // probeField();
    galaxySimulationP3M(int(5e4), 200);
    // galaxySimulationPM(int(5e4), 200);
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
  }
}

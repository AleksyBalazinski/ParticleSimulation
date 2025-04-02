#include <iostream>
#include "demos.h"

int main() {
  try {
    // smallSimPP();
    // bigSimPP();
    // smallSimP3M();
    // probeField();
    // galaxySimulationP3M(int(5e4), 200);
    galaxySimulationPM(int(5e4), 200);
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
  }
}

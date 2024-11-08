#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include "ppMethod.h"
#include "utils.h"

int main() {
  auto posRange = std::make_pair(Vec3(-100.0, -100.0, -100.0), Vec3(100.0, 100.0, 100.0));
  auto vRange = std::make_pair(Vec3(-20.0, -20.0, -20.0), Vec3(20.0, 20.0, 20.0));
  auto massRange = std::make_pair(5, 10);

  const int N = 3;
  auto state = randomInitialState(N, posRange, vRange);
  auto masses = randomMasses(N, massRange);

  std::string savePath = ppMethod(state, masses, 5.0);
  std::cout << "saving at " << savePath << '\n';
}

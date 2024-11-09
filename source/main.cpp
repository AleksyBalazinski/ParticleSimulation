#include <iostream>
#include <string>
#include "ppMethod.h"
#include "utils.h"

int main() {
  auto posRange = std::make_pair(Vec3(-10.0, -10.0, -10.0), Vec3(10.0, 10.0, 10.0));
  auto vRange = std::make_pair(Vec3(-2.0, -2.0, -2.0), Vec3(2.0, 2.0, 0.0));
  auto massRange = std::make_pair(20, 30);

  const int N = 3;
  auto state = randomInitialState(N, posRange, vRange);
  auto masses = randomMasses(N, massRange);

  std::string saveDirPath = ppMethod(state, masses, 10.0);
  std::cout << "saving at " << saveDirPath << '\n';
}

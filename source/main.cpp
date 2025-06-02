#include <chrono>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include "demos.h"
using namespace std;
using namespace std::chrono;
#ifdef CUDA
#include "gpu.h"
#endif

void printInfo(int argc, char* argv[]) {
  std::cout << "running " << argv[1] << " with arguments: ";
  for (int i = 2; i < argc; ++i) {
    std::cout << argv[i];
    if (i < argc - 1)
      std::cout << ", ";
  }
  std::cout << '\n';
}

void handleOptions(int argc, char* argv[]) {
  if (argc < 2) {
    throw std::runtime_error("Unknown option");
  }

  printInfo(argc, argv);

  std::string option(argv[1]);
  std::string arg2Str = (argc > 2) ? argv[2] : "";
  const char* arg2 = arg2Str.c_str();

  const std::unordered_map<std::string, std::function<void()>> optionHandlers = {
      {"anisotropy", [&]() { anisotropy(arg2); }},
      {"finite-diff-err", [&]() { finiteDiffRelError(arg2); }},
      {"assignment-schemes", [&]() { assignmentSchemes(arg2); }},
      {"poor-man-vs-laplacian", [&]() { poorLaplace(arg2); }},
      {"pm-accuracy-soft", [&]() { pmAccuracy(arg2, true); }},
      {"pm-accuracy-no-soft", [&]() { pmAccuracy(arg2, false); }},
      {"pm-optimal", [&]() { pmOptimal(arg2); }},
      {"pm-optimal-var-diameter", [&]() { pmOptimalVaryingDiameter(arg2); }},
      {"pm-optimal-assignments", [&]() { pmOptimalAssignments(arg2); }},
      {"p3m-accuracy-assignments", [&]() { p3mAccuracyAssignments(arg2); }},
      {"p3m-accuracy-shapes", [&]() { p3mAccuracyShapes(arg2); }},
      {"bh-accuracy", [&]() { bhAccuracy(arg2); }},
      {"galaxy-sim-pm", [&]() { galaxySimulationPM(arg2); }},
      {"galaxy-sim-p3m", [&]() { galaxySimulationP3M(arg2); }},
      {"galaxy-sim-bh", [&]() { galaxySimulationBH(arg2); }},
      {"cluster-sim-bh", [&]() { plummerBH(arg2); }},
      {"galaxy-collision-sim-bh", [&]() { galaxyCollisionBH(arg2); }}};

  auto it = optionHandlers.find(option);
  if (it != optionHandlers.end()) {
    it->second();
  } else {
    throw std::runtime_error("Unknown option: " + option);
  }
}

int main(int argc, char* argv[]) {
#ifdef CUDA
  printCudaVersion();
#endif
  try {
    handleOptions(argc, argv);
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
  }
}
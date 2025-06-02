#include "FFTWAdapter.h"
#ifdef CUDA
#include "PMMethodGPU.h"
#endif
#include <iostream>
#include "barnesHut.h"
#include "demos.h"
#include "diskSampler.h"
#include "diskSamplerLinear.h"
#include "externalFields.h"
#include "grid.h"
#include "kissFFTAdapter.h"
#include "p3mMethod.h"
#include "plummerSampler.h"
#include "pmMethod.h"
#include "ppMethod.h"
#include "utils.h"

std::vector<Vec3> linearSpacing(int n, float dx, Vec3 center, Vec3 dir) {
  std::vector<Vec3> state(2 * n);

  for (int i = 0; i < n; ++i) {
    state[i] = center + i * dx * dir;
    state[n + i] = Vec3::zero();
  }

  return state;
}

template <typename Container1, typename Container2>
auto cartesianProduct(const Container1& c1, const Container2& c2) {
  using T1 = typename Container1::value_type;
  using T2 = typename Container2::value_type;
  std::vector<std::pair<T1, T2>> result;

  for (const auto& a : c1) {
    for (const auto& b : c2) {
      result.emplace_back(a, b);
    }
  }

  return result;
}

template <typename Container1, typename Container2, typename Container3>
auto cartesianProduct(const Container1& c1, const Container2& c2, const Container3& c3) {
  using T1 = typename Container1::value_type;
  using T2 = typename Container2::value_type;
  using T3 = typename Container3::value_type;
  std::vector<std::tuple<T1, T2, T3>> result;

  for (const auto& a : c1) {
    for (const auto& b : c2) {
      for (const auto& c : c3) {
        result.emplace_back(a, b, c);
      }
    }
  }

  return result;
}

std::string toString(InterpolationScheme scheme) {
  switch (scheme) {
    case InterpolationScheme::NGP:
      return "NGP";
    case InterpolationScheme::CIC:
      return "CIC";
    case InterpolationScheme::TSC:
      return "TSC";
    default:
      return "Unknown";
  }
}

std::string toString(FiniteDiffScheme scheme) {
  switch (scheme) {
    case FiniteDiffScheme::TWO_POINT:
      return "2-point";
    case FiniteDiffScheme::FOUR_POINT:
      return "4-point";
    default:
      return "Unknown";
  }
}

std::string toString(GreensFunction gf) {
  switch (gf) {
    case GreensFunction::S1_OPTIMAL:
      return "S1";
    case GreensFunction::S2_OPTIMAL:
      return "S2";
    default:
      return "Unknown";
  }
}

CloudShape cloudFromGreensFunc(GreensFunction gf) {
  switch (gf) {
    case GreensFunction::S1_OPTIMAL:
      return CloudShape::S1;
    case GreensFunction::S2_OPTIMAL:
      return CloudShape::S2;
    default:
      throw std::exception("no corresponding particle shape");
  }
}

void anisotropy(const char* outputDir) {
  /*
  units:
  [M] = G(solar mass)
  [L] = kpc
  [T] = Myr
  */
  const int n = 400;
  const float G = 4.5e-3f;

  std::vector<float> masses(n, 0);
  masses[0] = 1.0f;

  auto externalField = [](Vec3 pos) -> Vec3 { return Vec3::zero(); };
  auto externalPotential = [](Vec3 pos) -> float { return 0; };

  auto gridPoints = std::make_tuple(64, 64, 32);
  std::array<int, 3> dims = {std::get<2>(gridPoints), std::get<1>(gridPoints),
                             std::get<0>(gridPoints)};
  FFTWAdapter fftAdapter(dims);
  Grid grid(gridPoints, fftAdapter);
  auto effectiveBoxSize = std::make_tuple(60.0f, 60.0f, 30.0f);
  float H = std::get<0>(effectiveBoxSize) / (std::get<0>(gridPoints) / 2);
  float DT = 1;
  float particleDiameter = 0;
  int simLength = 1;

  {
    StateRecorder stateRecorder(n, simLength + 1, outputDir, "positions.dat", "", "", "", "",
                                "field.dat");

    std::vector<Vec3> state =
        linearSpacing(n, 0.075f /*dx*/, Vec3(30, 30, 15) /*center*/, Vec3(1, 0, 0));

    PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
                InterpolationScheme::TSC, FiniteDiffScheme::TWO_POINT,
                GreensFunction::DISCRETE_LAPLACIAN, particleDiameter, grid);

    pm.run(stateRecorder, simLength, false /*diagnostics*/, true /*field*/);
  }

  {
    StateRecorder stateRecorder(n, simLength + 1, outputDir, "positions_offset.dat", "", "", "", "",
                                "field_offset.dat");

    std::vector<Vec3> state =
        linearSpacing(n, 0.075f /*dx*/, Vec3(30 - H / 2, 30, 15) /*center*/, Vec3(1, 0, 0));

    PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
                InterpolationScheme::TSC, FiniteDiffScheme::TWO_POINT,
                GreensFunction::DISCRETE_LAPLACIAN, particleDiameter, grid);

    pm.run(stateRecorder, simLength, false /*diagnostics*/, true /*field*/);
  }

  {
    StateRecorder stateRecorder(n, simLength + 1, outputDir, "positions_diag.dat", "", "", "", "",
                                "field_diag.dat");

    std::vector<Vec3> state = linearSpacing(n, 0.075f /*dx*/, Vec3(30, 30, 15) /*center*/,
                                            1 / std::sqrtf(2) * Vec3(1, 1, 0));

    PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
                InterpolationScheme::TSC, FiniteDiffScheme::TWO_POINT,
                GreensFunction::DISCRETE_LAPLACIAN, particleDiameter, grid);

    pm.run(stateRecorder, simLength, false /*diagnostics*/, true /*field*/);
  }
}

void finiteDiffRelError(const char* outputDir) {
  /*
  units:
  [M] = G(solar mass)
  [L] = kpc
  [T] = Myr
  */
  const int n = 400;
  const float G = 4.5e-3f;

  std::vector<float> masses(n, 0);
  masses[0] = 1.0f;

  auto externalField = [](Vec3 pos) -> Vec3 { return Vec3::zero(); };
  auto externalPotential = [](Vec3 pos) -> float { return 0; };

  auto gridPoints = std::make_tuple(64, 64, 32);
  std::array<int, 3> dims = {std::get<2>(gridPoints), std::get<1>(gridPoints),
                             std::get<0>(gridPoints)};
  FFTWAdapter fftAdapter(dims);
  Grid grid(gridPoints, fftAdapter);
  auto effectiveBoxSize = std::make_tuple(60.0f, 60.0f, 30.0f);
  float H = std::get<0>(effectiveBoxSize) / (std::get<0>(gridPoints) / 2);
  float DT = 1;
  float particleDiameter = 0;
  int simLength = 1;

  {
    StateRecorder stateRecorder(n, simLength + 1, outputDir, "positions-2.dat", "", "", "", "",
                                "field-2.dat");

    std::vector<Vec3> state =
        linearSpacing(n, 0.075f /*dx*/, Vec3(30, 30, 15) /*center*/, Vec3(1, 0, 0));

    PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
                InterpolationScheme::TSC, FiniteDiffScheme::TWO_POINT,
                GreensFunction::DISCRETE_LAPLACIAN, particleDiameter, grid);

    pm.run(stateRecorder, simLength, false /*diagnostics*/, true /*field*/);
  }

  {
    StateRecorder stateRecorder(n, simLength + 1, outputDir, "positions-4.dat", "", "", "", "",
                                "field-4.dat");

    std::vector<Vec3> state =
        linearSpacing(n, 0.075f /*dx*/, Vec3(30, 30, 15) /*center*/, Vec3(1, 0, 0));

    PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
                InterpolationScheme::TSC, FiniteDiffScheme::FOUR_POINT,
                GreensFunction::DISCRETE_LAPLACIAN, particleDiameter, grid);

    pm.run(stateRecorder, simLength, false /*diagnostics*/, true /*field*/);
  }
}

void assignmentSchemes(const char* outputDir) {
  /*
  units:
  [M] = G(solar mass)
  [L] = kpc
  [T] = Myr
  */
  const int n = 400;
  const float G = 4.5e-3f;

  std::vector<float> masses(n, 0);
  masses[0] = 1.0f;

  auto externalField = [](Vec3 pos) -> Vec3 { return Vec3::zero(); };
  auto externalPotential = [](Vec3 pos) -> float { return 0; };

  auto gridPoints = std::make_tuple(64, 64, 32);
  std::array<int, 3> dims = {std::get<2>(gridPoints), std::get<1>(gridPoints),
                             std::get<0>(gridPoints)};
  FFTWAdapter fftAdapter(dims);
  Grid grid(gridPoints, fftAdapter);
  auto effectiveBoxSize = std::make_tuple(60.0f, 60.0f, 30.0f);
  float H = std::get<0>(effectiveBoxSize) / (std::get<0>(gridPoints) / 2);
  float DT = 1;
  float particleDiameter = 0;
  int simLength = 1;

  {
    StateRecorder stateRecorder(n, simLength + 1, outputDir, "positions-tsc.dat", "", "", "", "",
                                "field-tsc.dat");

    std::vector<Vec3> state =
        linearSpacing(n, 0.075f /*dx*/, Vec3(30, 30, 15) /*center*/, Vec3(1, 0, 0));

    PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
                InterpolationScheme::TSC, FiniteDiffScheme::FOUR_POINT,
                GreensFunction::DISCRETE_LAPLACIAN, particleDiameter, grid);

    pm.run(stateRecorder, simLength, false /*diagnostics*/, true /*field*/);
  }

  {
    StateRecorder stateRecorder(n, simLength + 1, outputDir, "positions-cic.dat", "", "", "", "",
                                "field-cic.dat");

    std::vector<Vec3> state =
        linearSpacing(n, 0.075f /*dx*/, Vec3(30, 30, 15) /*center*/, Vec3(1, 0, 0));

    PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
                InterpolationScheme::CIC, FiniteDiffScheme::FOUR_POINT,
                GreensFunction::DISCRETE_LAPLACIAN, particleDiameter, grid);

    pm.run(stateRecorder, simLength, false /*diagnostics*/, true /*field*/);
  }

  {
    StateRecorder stateRecorder(n, simLength + 1, outputDir, "positions-ngp.dat", "", "", "", "",
                                "field-ngp.dat");

    std::vector<Vec3> state =
        linearSpacing(n, 0.075f /*dx*/, Vec3(30, 30, 15) /*center*/, Vec3(1, 0, 0));

    PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
                InterpolationScheme::NGP, FiniteDiffScheme::FOUR_POINT,
                GreensFunction::DISCRETE_LAPLACIAN, particleDiameter, grid);

    pm.run(stateRecorder, simLength, false /*diagnostics*/, true /*field*/);
  }
}

void poorLaplace(const char* outputDir) {
  /*
  units:
  [M] = G(solar mass)
  [L] = kpc
  [T] = Myr
  */
  const int n = 400;
  const float G = 4.5e-3f;

  std::vector<float> masses(n, 0);
  masses[0] = 1.0f;

  auto externalField = [](Vec3 pos) -> Vec3 { return Vec3::zero(); };
  auto externalPotential = [](Vec3 pos) -> float { return 0; };

  auto gridPoints = std::make_tuple(64, 64, 32);
  std::array<int, 3> dims = {std::get<2>(gridPoints), std::get<1>(gridPoints),
                             std::get<0>(gridPoints)};
  FFTWAdapter fftAdapter(dims);
  Grid grid(gridPoints, fftAdapter);
  auto effectiveBoxSize = std::make_tuple(60.0f, 60.0f, 30.0f);
  float H = std::get<0>(effectiveBoxSize) / (std::get<0>(gridPoints) / 2);
  float DT = 1;
  float particleDiameter = 0;
  int simLength = 1;

  {
    StateRecorder stateRecorder(n, simLength + 1, outputDir, "positions-poor-man.dat", "", "", "",
                                "", "field-poor-man.dat");

    std::vector<Vec3> state =
        linearSpacing(n, 0.075f /*dx*/, Vec3(30, 30, 15) /*center*/, Vec3(1, 0, 0));

    PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
                InterpolationScheme::TSC, FiniteDiffScheme::FOUR_POINT, GreensFunction::POOR_MAN,
                particleDiameter, grid);

    pm.run(stateRecorder, simLength, false /*diagnostics*/, true /*field*/);
  }

  {
    StateRecorder stateRecorder(n, simLength + 1, outputDir, "positions-discrete-lap.dat", "", "",
                                "", "", "field-discrete-lap.dat");

    std::vector<Vec3> state =
        linearSpacing(n, 0.075f /*dx*/, Vec3(30, 30, 15) /*center*/, Vec3(1, 0, 0));

    PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
                InterpolationScheme::TSC, FiniteDiffScheme::FOUR_POINT,
                GreensFunction::DISCRETE_LAPLACIAN, particleDiameter, grid);

    pm.run(stateRecorder, simLength, false /*diagnostics*/, true /*field*/);
  }
}

void pmAccuracy(const char* outputDir, bool ppSoftening) {
  Vec3 galaxyCenter = Vec3(30, 30, 15);
  float rb = 3.0f;
  float mb = 60.0f;
  float rd = 15.0f;
  float md = 15.0f;
  float thickness = 0.3f;
  float G = 4.5e-3f;
  int n = int(1e4);

  DiskSamplerLinear diskSampler(42);
  std::vector<Vec3> state = diskSampler.sample(galaxyCenter, rb, mb, rd, md, thickness, G, n);
  std::vector<float> masses(n, md / n);
  auto externalField = [](Vec3 pos) -> Vec3 { return Vec3::zero(); };
  auto externalPotential = [](Vec3 pos) -> float { return 0; };

  float DT = 1;
  int simLength = 0;
  auto effectiveBoxSize = std::make_tuple(60.0f, 60.0f, 30.0f);
  bool recordField = true;
  int NgStart = 30;
  int NgEnd = 70;
  int NgStep = 4;
  int NgPPStart = NgStart;
  int NgPPEnd = NgEnd;
  if (ppSoftening == false) {
    NgPPStart = 0;
    NgPPEnd = 0;
  }
  for (int Ng = NgPPStart; Ng <= NgPPEnd; Ng += NgStep) {
    float eps = ppSoftening ? std::get<0>(effectiveBoxSize) / Ng : 0.0f;
    std::string forceFile = "force-pp-" + std::to_string(Ng) + ".dat";
    StateRecorder stateRecorder(n, simLength + 1, outputDir, "positions-pp.dat", "", "", "", "",
                                forceFile.c_str());
    ppMethodLeapfrog(state, masses, simLength, DT, G, eps, stateRecorder, recordField);
  }

  std::array<InterpolationScheme, 3> interpolationSchemes = {
      InterpolationScheme::NGP, InterpolationScheme::CIC, InterpolationScheme::TSC};

  std::array<FiniteDiffScheme, 2> finiteDiffSchemes = {FiniteDiffScheme::TWO_POINT,
                                                       FiniteDiffScheme::FOUR_POINT};

  auto allCombinations = cartesianProduct(interpolationSchemes, finiteDiffSchemes);

  for (auto [is, fds] : allCombinations) {
    for (int Ng = NgStart; Ng <= NgEnd; Ng += NgStep) {
      auto gridPoints = std::make_tuple(Ng * 2, Ng * 2, Ng);
      std::array<int, 3> dims = {std::get<2>(gridPoints), std::get<1>(gridPoints),
                                 std::get<0>(gridPoints)};

      float H = std::get<0>(effectiveBoxSize) / (std::get<0>(gridPoints) / 2);
      FFTWAdapter fftAdapter(dims);
      Grid grid(gridPoints, fftAdapter);
      float particleDiameter = 0;
      std::string forceFile =
          "force-pm-" + toString(is) + "-" + toString(fds) + "-" + std::to_string(Ng) + ".dat";
      std::cout << forceFile << '\n';
      StateRecorder stateRecorder(n, simLength + 1, outputDir, "positions-pm.dat", "", "", "", "",
                                  forceFile.c_str());

      PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G, is,
                  fds, GreensFunction::DISCRETE_LAPLACIAN, particleDiameter, grid);

      pm.run(stateRecorder, simLength, false /*diagnostics*/, recordField);
    }
  }
}

void pmOptimal(const char* outputDir) {
  /*
  units:
  [M] = G(solar mass)
  [L] = kpc
  [T] = Myr
  */
  const int n = 400;
  const float G = 4.5e-3f;

  std::vector<float> masses(n, 0);
  masses[0] = 1.0f;

  auto externalField = [](Vec3 pos) -> Vec3 { return Vec3::zero(); };
  auto externalPotential = [](Vec3 pos) -> float { return 0; };

  auto gridPoints = std::make_tuple(64, 64, 32);
  std::array<int, 3> dims = {std::get<2>(gridPoints), std::get<1>(gridPoints),
                             std::get<0>(gridPoints)};
  FFTWAdapter fftAdapter(dims);
  Grid grid(gridPoints, fftAdapter);
  auto effectiveBoxSize = std::make_tuple(60.0f, 60.0f, 30.0f);
  float H = std::get<0>(effectiveBoxSize) / (std::get<0>(gridPoints) / 2);
  float DT = 1;
  float particleDiameter = 4 * H;
  int simLength = 1;

  {
    StateRecorder stateRecorder(n, simLength + 1, outputDir, "positions-s1.dat", "", "", "", "",
                                "field-s1.dat");

    std::vector<Vec3> state =
        linearSpacing(n, 0.075f /*dx*/, Vec3(30, 30, 15) /*center*/, Vec3(1, 0, 0));

    PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
                InterpolationScheme::TSC, FiniteDiffScheme::FOUR_POINT, GreensFunction::S1_OPTIMAL,
                particleDiameter, grid);

    pm.run(stateRecorder, simLength, false /*diagnostics*/, true /*field*/);
  }

  {
    StateRecorder stateRecorder(n, simLength + 1, outputDir, "positions-s2.dat", "", "", "", "",
                                "field-s2.dat");

    std::vector<Vec3> state =
        linearSpacing(n, 0.075f /*dx*/, Vec3(30, 30, 15) /*center*/, Vec3(1, 0, 0));

    PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
                InterpolationScheme::TSC, FiniteDiffScheme::FOUR_POINT, GreensFunction::S2_OPTIMAL,
                particleDiameter, grid);

    pm.run(stateRecorder, simLength, false /*diagnostics*/, true /*field*/);
  }
}

void pmOptimalVaryingDiameter(const char* outputDir) {
  /*
  units:
  [M] = G(solar mass)
  [L] = kpc
  [T] = Myr
  */
  const int n = 400;
  const float G = 4.5e-3f;

  std::vector<float> masses(n, 0);
  masses[0] = 1.0f;

  auto externalField = [](Vec3 pos) -> Vec3 { return Vec3::zero(); };
  auto externalPotential = [](Vec3 pos) -> float { return 0; };

  auto gridPoints = std::make_tuple(64, 64, 32);
  std::array<int, 3> dims = {std::get<2>(gridPoints), std::get<1>(gridPoints),
                             std::get<0>(gridPoints)};
  FFTWAdapter fftAdapter(dims);
  Grid grid(gridPoints, fftAdapter);
  auto effectiveBoxSize = std::make_tuple(60.0f, 60.0f, 30.0f);
  float H = std::get<0>(effectiveBoxSize) / (std::get<0>(gridPoints) / 2);
  float DT = 1;

  int simLength = 1;

  for (int i = 5; i <= 30; i++) {
    float particleDiameter = 0.2f * i * H;

    std::string fieldFilename = "field-s1-" + std::to_string(i) + ".dat";
    StateRecorder stateRecorder(n, simLength + 1, outputDir, "positions-s1.dat", "", "", "", "",
                                fieldFilename.c_str());

    std::vector<Vec3> state =
        linearSpacing(n, 0.075f /*dx*/, Vec3(30, 30, 15) /*center*/, Vec3(1, 0, 0));

    PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
                InterpolationScheme::TSC, FiniteDiffScheme::FOUR_POINT, GreensFunction::S1_OPTIMAL,
                particleDiameter, grid);

    pm.run(stateRecorder, simLength, false /*diagnostics*/, true /*field*/);
  }

  for (int i = 5; i <= 30; ++i) {
    float particleDiameter = 0.2f * i * H;

    std::string fieldFilename = "field-s2-" + std::to_string(i) + ".dat";
    StateRecorder stateRecorder(n, simLength + 1, outputDir, "positions-s2.dat", "", "", "", "",
                                fieldFilename.c_str());

    std::vector<Vec3> state =
        linearSpacing(n, 0.075f /*dx*/, Vec3(30, 30, 15) /*center*/, Vec3(1, 0, 0));

    PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
                InterpolationScheme::TSC, FiniteDiffScheme::FOUR_POINT, GreensFunction::S2_OPTIMAL,
                particleDiameter, grid);

    pm.run(stateRecorder, simLength, false /*diagnostics*/, true /*field*/);
  }
}

void pmOptimalAssignments(const char* outputDir) {
  /*
  units:
  [M] = G(solar mass)
  [L] = kpc
  [T] = Myr
  */
  const int n = 400;
  const float G = 4.5e-3f;

  std::vector<float> masses(n, 0);
  masses[0] = 1.0f;
  std::vector<Vec3> state =
      linearSpacing(n, 0.075f /*dx*/, Vec3(30, 30, 15) /*center*/, Vec3(1, 0, 0));

  auto externalField = [](Vec3 pos) -> Vec3 { return Vec3::zero(); };
  auto externalPotential = [](Vec3 pos) -> float { return 0; };

  auto gridPoints = std::make_tuple(64, 64, 32);
  std::array<int, 3> dims = {std::get<2>(gridPoints), std::get<1>(gridPoints),
                             std::get<0>(gridPoints)};
  FFTWAdapter fftAdapter(dims);
  Grid grid(gridPoints, fftAdapter);
  auto effectiveBoxSize = std::make_tuple(60.0f, 60.0f, 30.0f);
  float H = std::get<0>(effectiveBoxSize) / (std::get<0>(gridPoints) / 2);
  float DT = 1;
  float particleDiameter = 4 * H;
  int simLength = 1;

  std::array<InterpolationScheme, 3> interpolationSchemes = {
      InterpolationScheme::NGP, InterpolationScheme::CIC, InterpolationScheme::TSC};

  for (const auto& is : interpolationSchemes) {
    std::string fieldFile = "field-s1-" + toString(is) + ".dat";
    StateRecorder stateRecorder(n, simLength + 1, outputDir, "positions-s1.dat", "", "", "", "",
                                fieldFile.c_str());

    PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G, is,
                FiniteDiffScheme::FOUR_POINT, GreensFunction::S1_OPTIMAL, particleDiameter, grid);

    pm.run(stateRecorder, simLength, false /*diagnostics*/, true /*field*/);
  }
}

void p3mAccuracyAssignments(const char* outputDir) {
  Vec3 galaxyCenter = Vec3(30, 30, 15);
  float rb = 3.0f;
  float mb = 60.0f;
  float rd = 15.0f;
  float md = 15.0f;
  float thickness = 0.3f;
  float G = 4.5e-3f;
  int n = int(1e4);

  DiskSamplerLinear diskSampler(42);
  std::vector<Vec3> state = diskSampler.sample(galaxyCenter, rb, mb, rd, md, thickness, G, n);
  std::vector<float> masses(n, md / n);
  auto externalField = [](Vec3 pos) -> Vec3 { return Vec3::zero(); };
  auto externalPotential = [](Vec3 pos) -> float { return 0; };

  float DT = 1;
  int simLength = 0;
  auto effectiveBoxSize = std::make_tuple(60.0f, 60.0f, 30.0f);
  bool recordField = true;
  int aMulStart = 5;
  int aMulEnd = 25;
  int aMulStep = 1;

  float eps = 0.0f;
  std::string forceFile = "force-pp.dat";
  StateRecorder stateRecorder(n, simLength + 1, outputDir, "positions-pp.dat", "", "", "", "",
                              forceFile.c_str());
  ppMethodLeapfrog(state, masses, simLength, DT, G, eps, stateRecorder, recordField);

  std::array<InterpolationScheme, 3> interpolationSchemes = {
      InterpolationScheme::NGP, InterpolationScheme::CIC, InterpolationScheme::TSC};

  std::array<FiniteDiffScheme, 2> finiteDiffSchemes = {FiniteDiffScheme::TWO_POINT,
                                                       FiniteDiffScheme::FOUR_POINT};

  auto allCombinations = cartesianProduct(interpolationSchemes, finiteDiffSchemes);

  for (auto [is, fds] : allCombinations) {
    for (int aMul = aMulStart; aMul <= aMulEnd; aMul += aMulStep) {
      auto gridPoints = std::make_tuple(128, 128, 64);
      std::array<int, 3> dims = {std::get<2>(gridPoints), std::get<1>(gridPoints),
                                 std::get<0>(gridPoints)};

      float H = std::get<0>(effectiveBoxSize) / (std::get<0>(gridPoints) / 2);
      float a = 0.2f * aMul * H;
      FFTWAdapter fftAdapter(dims);
      Grid grid(gridPoints, fftAdapter);
      std::string forceFile =
          "force-p3m-" + toString(is) + "-" + toString(fds) + "-" + std::to_string(aMul) + ".dat";
      std::cout << forceFile << '\n';
      StateRecorder stateRecorder(n, simLength + 1, outputDir, "positions-pm.dat", "", "", "", "",
                                  forceFile.c_str());

      PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G, is,
                  fds, GreensFunction::S1_OPTIMAL, a, grid);

      float re = a;
      float softeningLength = 0.0f;
      P3MMethod p3m(pm, effectiveBoxSize, re, a, H, softeningLength, CloudShape::S1, false, false);

      p3m.run(stateRecorder, simLength, false, recordField);
    }
  }
}

void p3mAccuracyShapes(const char* outputDir) {
  Vec3 galaxyCenter = Vec3(30, 30, 15);
  float rb = 3.0f;
  float mb = 60.0f;
  float rd = 15.0f;
  float md = 15.0f;
  float thickness = 0.3f;
  float G = 4.5e-3f;
  int n = int(1e4);

  DiskSamplerLinear diskSampler(42);
  std::vector<Vec3> state = diskSampler.sample(galaxyCenter, rb, mb, rd, md, thickness, G, n);
  std::vector<float> masses(n, md / n);
  auto externalField = [](Vec3 pos) -> Vec3 { return Vec3::zero(); };
  auto externalPotential = [](Vec3 pos) -> float { return 0; };

  float DT = 1;
  int simLength = 0;
  auto effectiveBoxSize = std::make_tuple(60.0f, 60.0f, 30.0f);
  bool recordField = true;
  int aMulStart = 5;
  int aMulEnd = 25;
  int aMulStep = 1;

  float eps = 0.0f;
  std::string forceFile = "force-pp.dat";
  StateRecorder stateRecorder(n, simLength + 1, outputDir, "positions-pp.dat", "", "", "", "",
                              forceFile.c_str());
  ppMethodLeapfrog(state, masses, simLength, DT, G, eps, stateRecorder, recordField);

  std::array<FiniteDiffScheme, 2> finiteDiffSchemes = {FiniteDiffScheme::TWO_POINT,
                                                       FiniteDiffScheme::FOUR_POINT};

  std::array<GreensFunction, 2> greensFunctions = {GreensFunction::S1_OPTIMAL,
                                                   GreensFunction::S2_OPTIMAL};

  auto allCombinations = cartesianProduct(greensFunctions, finiteDiffSchemes);

  for (auto [gf, fds] : allCombinations) {
    for (int aMul = aMulStart; aMul <= aMulEnd; aMul += aMulStep) {
      auto gridPoints = std::make_tuple(128, 128, 64);
      std::array<int, 3> dims = {std::get<2>(gridPoints), std::get<1>(gridPoints),
                                 std::get<0>(gridPoints)};

      float H = std::get<0>(effectiveBoxSize) / (std::get<0>(gridPoints) / 2);
      float a = 0.2f * aMul * H;
      FFTWAdapter fftAdapter(dims);
      Grid grid(gridPoints, fftAdapter);
      std::string forceFile = "force-p3m-TSC-" + toString(fds) + "-" + toString(gf) + "-" +
                              std::to_string(aMul) + ".dat";
      std::cout << forceFile << '\n';
      StateRecorder stateRecorder(n, simLength + 1, outputDir, "positions-pm.dat", "", "", "", "",
                                  forceFile.c_str());

      PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
                  InterpolationScheme::TSC, fds, gf, a, grid);

      float re = a;
      float softeningLength = 0.0f;
      P3MMethod p3m(pm, effectiveBoxSize, re, a, H, softeningLength, cloudFromGreensFunc(gf), false,
                    false);

      p3m.run(stateRecorder, simLength, false, recordField);
    }
  }
}

void galaxySimulationPM(const char* outputDir) {
  /*
    units:
    [M] = G(solar mass)
    [L] = kpc
    [T] = Myr
    */
  Vec3 galaxyCenter = Vec3(30, 30, 15);
  float rb = 3.0f;
  float mb = 60.0f;
  float rd = 15.0f;
  float md = 15.0f;
  float thickness = 0.3f;
  float G = 4.5e-3f;
  int n = int(5e4);
  int simLength = 200;

  DiskSamplerLinear diskSampler(42);
  std::vector<Vec3> state = diskSampler.sample(galaxyCenter, rb, mb, rd, md, thickness, G, n);

  std::vector<float> masses(n, md / n);
  auto externalField = [galaxyCenter, rb, mb, G](Vec3 pos) -> Vec3 {
    return sphRadDecrField(pos, galaxyCenter, rb, mb, G);
  };
  auto externalPotential = [galaxyCenter, rb, mb, G](Vec3 pos) -> float {
    return sphRadDecrFieldPotential(pos, galaxyCenter, rb, mb, G);
  };

  auto gridPoints = std::make_tuple(128, 128, 64);
  std::array<int, 3> dims = {std::get<2>(gridPoints), std::get<1>(gridPoints),
                             std::get<0>(gridPoints)};
  auto effectiveBoxSize = std::make_tuple(60.0f, 60.0f, 30.0f);
  float H = std::get<0>(effectiveBoxSize) / (std::get<0>(gridPoints) / 2);
  float DT = 1;

#ifdef CUDA
  PMMethodGPU pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
                 InterpolationScheme::TSC, FiniteDiffScheme::TWO_POINT,
                 GreensFunction::DISCRETE_LAPLACIAN, 0, gridPoints);
#else
  FFTWAdapter fftAdapter(dims);
  Grid grid(gridPoints, fftAdapter);
  PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
              InterpolationScheme::TSC, FiniteDiffScheme::TWO_POINT,
              GreensFunction::DISCRETE_LAPLACIAN, 0, grid);
#endif
  StateRecorder stateRecorder(n, simLength + 1, outputDir);
  pm.run(stateRecorder, simLength, true /*diagnostics*/);
}

void galaxySimulationP3M(const char* outputDir) {
  /*
    units:
    [M] = G(solar mass)
    [L] = kpc
    [T] = Myr
    */
  Vec3 galaxyCenter = Vec3(30, 30, 15);
  float rb = 3.0f;
  float mb = 60.0f;
  float rd = 15.0f;
  float md = 15.0f;
  float thickness = 0.3f;
  float G = 4.5e-3f;
  int simLength = 200;
  int n = int(5e4);

  DiskSamplerLinear diskSampler(42);
  std::vector<Vec3> state = diskSampler.sample(galaxyCenter, rb, mb, rd, md, thickness, G, n);

  std::vector<float> masses(n, md / n);
  auto externalField = [galaxyCenter, rb, mb, G](Vec3 pos) -> Vec3 {
    return sphRadDecrField(pos, galaxyCenter, rb, mb, G);
  };
  auto externalPotential = [galaxyCenter, rb, mb, G](Vec3 pos) -> float {
    return sphRadDecrFieldPotential(pos, galaxyCenter, rb, mb, G);
  };

  auto gridPoints = std::make_tuple(128, 128, 64);
  std::array<int, 3> dims = {std::get<2>(gridPoints), std::get<1>(gridPoints),
                             std::get<0>(gridPoints)};
  FFTWAdapter fftAdapter(dims);
  Grid grid(gridPoints, fftAdapter);
  auto effectiveBoxSize = std::make_tuple(60.0f, 60.0f, 30.0f);
  float H = std::get<0>(effectiveBoxSize) / (std::get<0>(gridPoints) / 2);
  float DT = 1;
  float a = 3 * H;

#ifdef CUDA
  PMMethodGPU pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
                 InterpolationScheme::TSC, FiniteDiffScheme::TWO_POINT,
                 GreensFunction::DISCRETE_LAPLACIAN, 0, gridPoints);
#else
  PMMethod pm(state, masses, effectiveBoxSize, externalField, externalPotential, H, DT, G,
              InterpolationScheme::TSC, FiniteDiffScheme::TWO_POINT, GreensFunction::S1_OPTIMAL, a,
              grid);
#endif

  float re = 0.7f * a;
  float softeningLength = 1.5f;
  P3MMethod p3m(pm, effectiveBoxSize, re, a, H, softeningLength, CloudShape::S1);

  StateRecorder stateRecorder(n, simLength + 1, outputDir);
  p3m.run(stateRecorder, simLength, true /*diagnostics*/);
}

void galaxySimulationBH(const char* outputDir) {
  Vec3 galaxyCenter = Vec3(30, 30, 15);
  float rb = 3.0f;
  float mb = 60.0f;
  float rd = 15.0f;
  float md = 15.0f;
  float thickness = 0.3f;
  float G = 4.5e-3f;
  int n = int(5e4);
  int simLength = 200;

  DiskSamplerLinear diskSampler(42);
  std::vector<Vec3> state = diskSampler.sample(galaxyCenter, rb, mb, rd, md, thickness, G, n);
  std::vector<float> masses(n, md / n);
  auto externalField = [galaxyCenter, rb, mb, G](Vec3 pos) -> Vec3 {
    return sphRadDecrField(pos, galaxyCenter, rb, mb, G);
  };
  auto externalPotential = [galaxyCenter, rb, mb, G](Vec3 pos) -> float {
    return sphRadDecrFieldPotential(pos, galaxyCenter, rb, mb, G);
  };

  Vec3 low(0, 0, 0);
  float H = 60;
  float DT = 1;
  float softeningLength = 1.5f;
  float theta = 1.0f;

  BH::BarnesHut bhSimulation(state, masses, externalField, externalPotential, low, H, G,
                             softeningLength, theta);

  StateRecorder stateRecorder(n, simLength + 1, outputDir);
  bhSimulation.run(stateRecorder, simLength, DT, true /*diagnostics*/);
}

void bhAccuracy(const char* outputDir) {
  Vec3 galaxyCenter = Vec3(30, 30, 15);
  float rb = 3.0f;
  float mb = 60.0f;
  float rd = 15.0f;
  float md = 15.0f;
  float thickness = 0.3f;
  float G = 4.5e-3f;
  int n = int(1e4);

  DiskSamplerLinear diskSampler(42);
  std::vector<Vec3> state = diskSampler.sample(galaxyCenter, rb, mb, rd, md, thickness, G, n);
  std::vector<float> masses(n, md / n);
  auto externalField = [](Vec3 pos) -> Vec3 { return Vec3::zero(); };
  auto externalPotential = [](Vec3 pos) -> float { return 0; };

  Vec3 low(0, 0, 0);
  float H = 60;
  float DT = 1;
  int simLength = 0;
  bool recordField = true;

  float softeningLength = 0.0f;
  std::string forceFile = "force-pp.dat";
  StateRecorder stateRecorder(n, simLength + 1, outputDir, "positions-pp.dat", "", "", "", "",
                              forceFile.c_str());
  ppMethodLeapfrog(state, masses, simLength, DT, G, softeningLength, stateRecorder, recordField);

  bool quadrupoles[] = {false, true};
  for (bool quadrupole : quadrupoles) {
    for (int thetaMul = 0; thetaMul < 12; thetaMul++) {
      float theta = 0.1f * thetaMul;
      std::string forceFile =
          "force-bh-" + std::to_string(thetaMul) + (quadrupole ? "-q" : "-m") + ".dat";
      StateRecorder stateRecorder(n, simLength + 1, outputDir, "positions-pm.dat", "", "", "", "",
                                  forceFile.c_str());

      BH::BarnesHut bhSimulation(state, masses, externalField, externalPotential, low, H, G,
                                 softeningLength, theta, quadrupole);

      bhSimulation.run(stateRecorder, simLength, DT, false, recordField);
    }
  }
}

void galaxyCollisionBH(const char* outputDir) {
  Vec3 galaxyCenter1 = Vec3(40, 30, 15);
  Vec3 galaxyCenter2 = Vec3(80, 30, 15);
  float rb = 3.0f;
  float mb = 60.0f;
  float rd = 15.0f;
  float md = 15.0f;
  float thickness = 0.3f;
  float G = 4.5e-3f;
  int n = int(3e4);

  DiskSamplerLinear diskSampler(42);
  std::vector<Vec3> state1 = diskSampler.sample(galaxyCenter1, rb, mb, rd, md, thickness, G, n, rb);
  std::vector<float> masses1(n, md / n);
  state1[0] = galaxyCenter1;
  state1[n] = Vec3::zero();
  masses1[0] = 60.0f;

  std::vector<Vec3> state2 = diskSampler.sample(galaxyCenter2, rb / 1.5f, mb / 1.5f, rd / 1.5f,
                                                md / 1.5f, thickness, G, n, rb / 1.5f);
  std::vector<float> masses2(n, md / n);
  state2[0] = galaxyCenter2;
  state2[n] = Vec3::zero();
  masses2[0] = 60.0f / 1.5f;

  std::vector<Vec3> state(4 * n);
  std::vector<float> masses(2 * n);
  for (int i = 0; i < n; ++i) {
    state[i] = state1[i];
    state[i + 2 * n] = state1[i + n];
    masses[i] = masses1[i];
  }
  for (int i = 0; i < n; ++i) {
    state[i + n] = state2[i];
    state[i + 3 * n] = state2[i + n];
    masses[i + n] = masses2[i];
  }

  auto externalField = [](Vec3 pos) -> Vec3 { return Vec3::zero(); };
  auto externalPotential = [](Vec3 pos) -> float { return 0; };

  Vec3 low(0, 0, 0);
  float H = 120;
  float DT = 1;
  int simLength = 200;
  float softeningLength = 1.3f;
  float theta = 1.0f;

  BH::BarnesHut bhSimulation(state, masses, externalField, externalPotential, low, H, G,
                             softeningLength, theta);

  StateRecorder stateRecorder(2 * n, simLength + 1, outputDir);
  bhSimulation.run(stateRecorder, simLength, DT, true /*diagnostics*/);
}

void plummerBH(const char* outputDir) {
  Vec3 clusterCenter = Vec3(30, 30, 30);
  float a = 2;
  float rMax = 15;
  float M = 1.0f;
  float G = 4.5e-3f;
  int n = int(1e4);

  PlummerSampler sampler(42);
  std::vector<Vec3> state = sampler.sample(clusterCenter, a, rMax, M, G, n);
  std::vector<float> masses(n, M / n);
  auto externalField = [](Vec3 pos) -> Vec3 { return Vec3::zero(); };
  auto externalPotential = [](Vec3 pos) -> float { return 0; };

  Vec3 low(0, 0, 0);
  float H = 60;
  float DT = 1;
  int simLength = 200;
  float softeningLength = 0.5f;
  float theta = 0.5f;

  BH::BarnesHut bhSimulation(state, masses, externalField, externalPotential, low, H, G,
                             softeningLength, theta);

  StateRecorder stateRecorder(n, simLength + 1, outputDir);
  bhSimulation.run(stateRecorder, simLength, DT, true /*diagnostics*/);
}
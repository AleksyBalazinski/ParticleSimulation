#include "common.h"
#include "disk_sampler_linear.cuh"
#include "pm.cuh"
#include "state_recorder.cuh"
#include "utils.cuh"
#include "vec3.cuh"

int main() {
  int simLength = 200;
  Vec3 galaxyCenter(30, 30, 30);
  float rb = 3.0f;
  float mb = 15.0f;
  float rd = 15.0f;
  float md = 60.0f;
  float thickness = 0.3f;
  float effectiveBoxSize = 60;
  float G = 4.5e-3f;

  float H = effectiveBoxSize / (NG / 2);
  float DT = 1;

  DiskSamplerLinear sampler;
  std::vector<Vec3> state = sampler.sample(galaxyCenter, rb, mb, rd, md, thickness, G, N);
  std::vector<float> masses(N, md / N);

  pmMethod(state, masses, effectiveBoxSize, H, DT, G, simLength);
}
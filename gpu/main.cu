#include "disk_sampler_linear.cuh"
#include "external_fields.cuh"
#include "pm.cuh"
#include "settings.cuh"
#include "state_recorder.cuh"
#include "triple.cuh"
#include "vec3.cuh"

int main() {
  int simLength = 200;
  Vec3 galaxyCenter(30, 30, 15);
  float rb = 3.0f;
  float mb = 60.0f;
  float rd = 15.0f;
  float md = 15.0f;
  float thickness = 0.3f;
  Triple<float> effectiveBoxSize(60, 60, 30);
  float G = 4.5e-3f;

  float H = effectiveBoxSize.x / (NGx / 2);
  float DT = 1;

  DiskSamplerLinear sampler;
  SphRadDecrFieldParams bulgeParams(galaxyCenter, rb, mb);
  std::vector<Vec3> state = sampler.sample(bulgeParams, rd, md, thickness, G, N);
  std::vector<float> masses(N, md / N);

  pmMethod(state, masses, effectiveBoxSize, H, DT, bulgeParams, G, simLength);
}
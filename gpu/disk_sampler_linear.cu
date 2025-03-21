#include "disk_sampler_linear.cuh"
#include "external_fields.cuh"
#include "utils.cuh"

Vec3 DiskSamplerLinear::getVelocity(Vec3 pos,
                                    SphRadDecrFieldParams bulgeParams,
                                    float rd,
                                    float md,
                                    float G) {
  Vec3 rVec = pos - bulgeParams.center;
  float r = rVec.getMagnitude();

  Vec3 rVecxy = rVec;
  rVecxy.z = 0;
  float rho = rVecxy.getMagnitude();

  float ra = r / rd;
  float sigma0 = 3 * md / (PI * rd * rd);
  float k = 2.5;
  float h = 0.66;
  float a = -k / (h * h);
  float gdVal = -G * sigma0 * (a * (ra - h) * (ra - h) + k);

  Vec3 gb = sphRadDecrFieldHost(pos, bulgeParams, G);
  Vec3 gd = gdVal * rVecxy / rVecxy.getMagnitude();
  Vec3 g = gb + gd;
  float gVal = g.getMagnitude();

  float v = std::sqrt(gVal * rho * rho / r);
  return Vec3(-v * rVec.y / rho, v * rVec.x / rho, 0);
}

float DiskSamplerLinear::sampleFromLinear(float rd) {
  std::uniform_real_distribution<float> u(0, 1);

  float cdf = u(re);
  auto f = [rd, cdf, this](float r) { return implicitInvCDFLinear(r, rd, cdf); };
  auto df = [rd, this](float r) { return implicitInvCDFLinearDer(r, rd); };
  return newton(f, df, rd / 2.0, 100, 0.01);
}

DiskSamplerLinear::DiskSamplerLinear() : re(std::random_device{}()) {}

std::vector<Vec3> DiskSamplerLinear::sample(SphRadDecrFieldParams bulgeParams,
                                            float rd,
                                            float md,
                                            float thickness,
                                            float G,
                                            int n) {
  std::uniform_real_distribution<float> u(0, 1);
  std::vector<Vec3> state(2 * n);

  for (int i = 0; i < n; ++i) {
    float phi = u(re) * 2 * PI;
    float r = sampleFromLinear(rd);

    float x = r * std::cos(phi);
    float y = r * std::sin(phi);
    float z = (thickness / 2) * (2 * u(re) - 1);

    state[i] = bulgeParams.center + Vec3(x, y, z);
    state[n + i] = getVelocity(state[i], bulgeParams, rd, md, G);
  }

  return state;
}

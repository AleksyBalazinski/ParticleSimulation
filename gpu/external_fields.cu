#include "external_fields.cuh"

__host__ __device__ SphRadDecrFieldParams::SphRadDecrFieldParams(Vec3 center, float r, float m)
    : center(center), r(r), m(m) {}

__device__ Vec3 sphRadDecrField(Vec3 pos, SphRadDecrFieldParams params, float G) {
  auto bulge = params.center;
  auto mb = params.m;
  auto rb = params.r;

  Vec3 d(pos.x - bulge.x, pos.y - bulge.y, pos.z - bulge.z);
  float r = norm3df(d.x, d.y, d.z);
  Vec3 dir(d.x / r, d.y / r, d.z / r);
  float g;
  if (r > rb) {
    g = -G * mb / (r * r);
  } else {
    g = -(G * mb / (rb * rb * rb)) * r * (4 - 3 * r / rb);
  }

  return Vec3(g * dir.x, g * dir.y, g * dir.z);
}

__host__ Vec3 sphRadDecrFieldHost(Vec3 pos, SphRadDecrFieldParams params, float G) {
  auto bulge = params.center;
  auto mb = params.m;
  auto rb = params.r;

  float r = (pos - bulge).getMagnitude();
  Vec3 dir = (pos - bulge) / r;
  float g;
  if (r > rb) {
    g = -G * mb / (r * r);
  } else {
    g = -(G * mb / std::powf(rb, 3)) * r * (4 - 3 * r / rb);
  }

  return g * dir;
}
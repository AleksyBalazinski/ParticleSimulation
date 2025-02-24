#include "ext_field_bind.cuh"

ExtFieldBind::ExtFieldBind(Vec3 center,
                           float rb,
                           float mb,
                           float G,
                           Vec3 (*field)(Vec3, Vec3, float, float, float))
    : center(center), rb(rb), mb(mb), G(G), field(field) {}

__host__ __device__ Vec3 ExtFieldBind::operator()(Vec3 pos) {
  return field(pos, this->center, this->rb, this->mb, this->G);
}
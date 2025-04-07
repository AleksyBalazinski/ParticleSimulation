template <typename T>
struct Triple {
  T x;
  T y;
  T z;
  __host__ __device__ Triple(T x, T y, T z) : x(x), y(y), z(z) {}
};

template <typename T>
__host__ __device__ Triple<T> makeTriple(T x, T y, T z) {
  return Triple<T>(x, y, z);
}
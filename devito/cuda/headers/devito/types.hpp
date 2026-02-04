#ifndef _DEVITO_CUDA_TYPES_H
#define _DEVITO_CUDA_TYPES_H

/* ============================================================================
 * Additional types and type helpers for the CUDA backend
 */
#include <cuda_runtime.h>
#include <utility>

namespace devito {
namespace cuda {

// Dynamic access to dim3 members
static inline int dim3_get(const dim3 &d, int rank) {
  if (rank == 0)
    return d.x;
  else if (rank == 1)
    return d.y;
  else if (rank == 2)
    return d.z;

  assert(false);
  return -1;
}

static void dim3_set(dim3 &d, int rank, int value) {
  switch (rank) {
  case 0:
    d.x = value;
    return;
  case 1:
    d.y = value;
    return;
  case 2:
    d.z = value;
    return;
  default:
    assert(false);
  }
}

struct dim4 {
  unsigned int w, x, y, z;
  __host__ __device__ dim4(const unsigned int vw = 1, const unsigned int vx = 1,
                           const unsigned int vy = 1, const unsigned int vz = 1)
      : w(vw), x(vx), y(vy), z(vz) {}
  __host__ __device__ dim4(uint3 v) : w(1), x(v.x), y(v.y), z(v.z) {}
  __host__ __device__ dim4(dim3 v) : w(1), x(v.x), y(v.y), z(v.z) {}
};

struct int3 {
  int x, y, z;
  __host__ __device__ int3(const int vx = 1, const int vy = 1, const int vz = 1)
      : x(vx), y(vy), z(vz) {}
  __host__ __device__ int3(uint3 v) : x(v.x), y(v.y), z(v.z) {}
  __host__ __device__ int3(dim3 v) : x(v.x), y(v.y), z(v.z) {}
};

struct int4 {
  int w, x, y, z;
  __host__ __device__ int4(const int vw = 1, const int vx = 1, const int vy = 1,
                           const int vz = 1)
      : w(vw), x(vx), y(vy), z(vz) {}
  __host__ __device__ int4(uint3 v) : w(1), x(v.x), y(v.y), z(v.z) {}
  __host__ __device__ int4(dim3 v) : w(1), x(v.x), y(v.y), z(v.z) {}
};

namespace pair_iterators {
template <typename T1, typename T2>
std::pair<T1, T2> operator++(std::pair<T1, T2> &it) {
  ++it.first;
  ++it.second;
  return it;
}
} // namespace pair_iterators

} // namespace cuda
} // namespace devito
#endif

#ifndef _DEVITO_CUDA_TYPES_H
#define _DEVITO_CUDA_TYPES_H

/* ============================================================================
 * Additional types and type helpers for the CUDA backend
 */
#include <cuda_runtime.h>

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
  __host__ __device__ dim4(unsigned int vw = 1, unsigned int vx = 1,
                           unsigned int vy = 1, unsigned int vz = 1)
      : w(vw), x(vx), y(vy), z(vz) {}
  __host__ __device__ dim4(uint3 v) : w(1), x(v.x), y(v.y), z(v.z) {}
  __host__ __device__ dim4(dim3 v) : w(1), x(v.x), y(v.y), z(v.z) {}
  __host__ __device__ operator uint3(void) const { return uint3{x, y, z}; }
};

} // namespace cuda
} // namespace devito
#endif

#ifndef _DEVITO_CUDA_UTIL_H
#define _DEVITO_CUDA_UTIL_H

#include <cassert>
#include <cuda_runtime.h>

namespace devito {
namespace cuda {

inline uint64_t next_pow2(uint64_t x) {
  return (__builtin_popcount(x) == 1 || x == 1)
             ? x
             : 1 << (64 - __builtin_clzl(x - 1));
}

static inline int ceil_div(int x, int y) { return (x + y - 1) / y; }
static inline int round_down_to_multiple(int x, int y) { return (x / y) * y; }

template <typename T> inline bool _cudaPointerIsAccessible(T *ptr) {
  struct cudaPointerAttributes attr = cudaPointerAttributes{};
  CudaChecked(cudaPointerGetAttributes(&attr, (const void *)ptr));
  return attr.devicePointer != nullptr;
}

template <typename T> inline bool _cudaPtrIsManaged(T *ptr) {
  struct cudaPointerAttributes attr = cudaPointerAttributes{};
  CudaChecked(cudaPointerGetAttributes(&attr, (const void *)ptr));
  return attr.type == cudaMemoryType::cudaMemoryTypeManaged;
}

template <typename T> inline bool _cudaPtrIsDeviceAccessible(T *ptr) {
  if (ptr == nullptr)
    return false;

  struct cudaPointerAttributes attr = cudaPointerAttributes{};
  auto ret = cudaPointerGetAttributes(&attr, (const void *)ptr);
  return ret == cudaSuccess && attr.devicePointer != nullptr;
}

template <typename T> inline bool _cudaPtrIsDevicePtr(T *ptr) {
  if (ptr == nullptr)
    return false;

  struct cudaPointerAttributes attr = cudaPointerAttributes{};
  auto ret = cudaPointerGetAttributes(&attr, (const void *)ptr);
  return ret == cudaSuccess &&
         attr.type == cudaMemoryType::cudaMemoryTypeDevice;
}

#define _setCudaConstant(TYPE, NAME, VAL)                                      \
  {                                                                            \
    TYPE tmp = VAL;                                                            \
    CudaChecked(cudaMemcpyToSymbol(NAME, &tmp, sizeof(TYPE)));                 \
  }

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

} // namespace cuda
} // namespace devito

#endif // _DEVITO_CUDA_UTIL_H

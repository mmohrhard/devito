#ifndef _DEVITO_CUDA_UTIL_H
#define _DEVITO_CUDA_UTIL_H

#include <cassert>
#include <cuda_runtime.h>
#include <devito/errors.hpp>

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

#define SFINAE_FOR_MEMBER_CHECK(NAME, MEMBER)                                      \
  template <typename T, typename = void> struct has_##NAME : std::false_type {}; \
  template <typename T>                                                        \
  struct has_##NAME<T, std::void_t<decltype(T::MEMBER)>>                    \
      : std::true_type {};                                                    \
  template <typename T> constexpr bool has_##NAME##_v = has_##NAME<T>::value;

// SFINAE to check for presence of host pointer member
// (Devito Arrays and dataobjs are nearly-interchangeable
// structs, but arrays don't have host pointers)
SFINAE_FOR_MEMBER_CHECK(host_ptr, data)

template <typename T>
typename std::enable_if<has_host_ptr_v<T>, void *>::type
select_ptr(T *obj, bool device) {
  if (device) {
    return obj->device_data;
  } else {
    return obj->data;
  }
}

template <typename T>
typename std::enable_if<!has_host_ptr_v<T>, void *>::type
select_ptr(T *obj, bool device) {
  if (device) {
    return obj->device_data;
  } else {
    // This should never happen, unless codegen has messed up
    assert(false && "Object does not have host_data member");
  }
}

// Upcoming savebuffer compression options
SFINAE_FOR_MEMBER_CHECK(is_compressed, compressed);
SFINAE_FOR_MEMBER_CHECK(compression_type, compression_type);
SFINAE_FOR_MEMBER_CHECK(compression_buf, compression_buf);

} // namespace cuda
} // namespace devito

#endif // _DEVITO_CUDA_UTIL_H

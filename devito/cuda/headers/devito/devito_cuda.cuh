#pragma once

#include "nvtx3/nvToolsExt.h"
#include <cmath>
#include <cuda.h>
#include <devito/jitify.hpp>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <nccl.h>
#include <set>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

#ifndef OPERATOR_STANDALONE
#include <Python.h>
#endif

#ifndef STRINGIFY
#define STRINGIFY(x) _stringify(x)
#define _stringify(x) #x
#endif

#ifndef NVRTC_CUDA_ARCH
// Default to compute capability 8.0 aka A100
#define NVRTC_CUDA_ARCH compute_80
#endif

#ifdef KERNEL_DEBUGGING
#define DEBUG_OPTS "-G",
#else
#define DEBUG_OPTS
#endif

#define NVRTC_OPTS                                                             \
  {"--ftz=true",                                                               \
   "--fmad=true",                                                              \
   "--prec-sqrt=true",                                                         \
   "--prec-div=true",                                                          \
   DEBUG_OPTS "--gpu-architecture=" STRINGIFY(NVRTC_CUDA_ARCH),                \
   "--std=c++11"}

/**
 * LOGGING
 */

enum LogLevel {
  CRITICAL = 50,
  FATAL = 50,
  ERROR = 40,
  WARNING = 30,
  WARN = 30,
  INFO = 20,
  DEBUG = 10,
  NOTSET = 0
};

typedef void (*LogHandler)(int logLevel, const char *message);

LogHandler _logHandler = nullptr;

extern "C" void setLogHandler(void *handler) {
  _logHandler = (LogHandler)handler;
}

template <typename... Args>
std::string string_format(const std::string &format, Args... args) {
  int size_s = std::snprintf(nullptr, 0, format.c_str(), args...) +
               1; // Extra space for '\0'
  if (size_s <= 0) {
    throw std::runtime_error("Error during formatting.");
  }
  auto size = static_cast<size_t>(size_s);
  std::unique_ptr<char[]> buf(new char[size]);
  std::snprintf(buf.get(), size, format.c_str(), args...);
  return std::string(buf.get(),
                     buf.get() + size - 1); // We don't want the '\0' inside
}

template <typename... Args>
inline void log(int logLevel, const std::string &format, Args... args) {
  std::string message = string_format(format, std::forward<Args>(args)...);

  if (_logHandler == nullptr)
    return;
  _logHandler(logLevel, message.c_str());
}

template <typename... Args>
inline void debug(const std::string &format, Args... args) {
  log(LogLevel::DEBUG, format, std::forward<Args>(args)...);
}

template <typename... Args>
inline void info(const std::string &format, Args... args) {
  log(LogLevel::INFO, format, std::forward<Args>(args)...);
}

template <typename... Args>
inline void warn(const std::string &format, Args... args) {
  log(LogLevel::WARNING, format, std::forward<Args>(args)...);
}

template <typename... Args>
inline void critical(const std::string &format, Args... args) {
  log(LogLevel::CRITICAL, format, std::forward<Args>(args)...);
}

#define Checked(...)                                                           \
  {                                                                            \
    auto _ret = (__VA_ARGS__);                                                 \
    if (_ret != 0) {                                                           \
      return _ret;                                                             \
    }                                                                          \
  }
/**
 * CUDA error helpers
 */
#define CudaChecked(f)                                                         \
  {                                                                            \
    if (!_cudaChecked((cudaError_t)(f), __FILE__, __LINE__))                   \
      return -1;                                                               \
  }
#define CudaCheckedEx(f, msg)                                                  \
  {                                                                            \
    if (!_cudaChecked((cudaError_t)(f), __FILE__, __LINE__, msg))              \
      return -1;                                                               \
  }
#define NcclChecked(f)                                                         \
  {                                                                            \
    if (!_ncclChecked((ncclResult_t)(f), __FILE__, __LINE__))                  \
      return -1;                                                               \
  }
#define NcclCheckedEx(f, msg)                                                  \
  {                                                                            \
    if (!_ncclChecked((ncclResult_t)(f), __FILE__, __LINE__, msg))             \
      return -1;                                                               \
  }

/**
 * Thread-local temporary variable helpers
 */
#define MAX_CUDA_DEVICES 16
#define PER_DEVICE_TEMP_GET(T, NAME, SIZE)                                     \
  T *NAME = nullptr;                                                           \
  static T *NAME##_device[MAX_CUDA_DEVICES] = {0};                             \
  {                                                                            \
    int device = 0;                                                            \
    CudaChecked(cudaGetDevice(&device));                                       \
    assert(device >= 0 && device < MAX_CUDA_DEVICES);                          \
    if (NAME##_device[device] == nullptr) {                                    \
      debug("allocating %llu bytes for " STRINGIFY(NAME) " on device %d",      \
            SIZE, device);                                                     \
      CudaChecked(cudaMallocAsync((void **)&NAME##_device[device], (SIZE),     \
                                  kernel_stream));                             \
      CudaChecked(                                                             \
          cudaMemsetAsync(NAME##_device[device], 1, (SIZE), kernel_stream));   \
    }                                                                          \
    NAME = NAME##_device[device];                                              \
  }

#define PER_DEVICE_TEMP_DESTROY(NAME)                                          \
  {                                                                            \
    int device = 0;                                                            \
    CudaChecked(cudaGetDevice(&device));                                       \
    assert(device >= 0 && device < MAX_CUDA_DEVICES);                          \
    CudaChecked(cudaFreeAsync(NAME##_device[device], kernel_stream));          \
    NAME##_device[device] = nullptr;                                           \
  }

#define PER_DEVICE_ARRAY_TEMP_DECLARE(NAME, ARRAYTYPE)                         \
  static ARRAYTYPE *NAME##_device[MAX_CUDA_DEVICES] = {0};                     \
  ARRAYTYPE *NAME = {0};

#define PER_DEVICE_ARRAY_TEMP_GET(NAME, DTYPE, ...)                            \
  {                                                                            \
    if (_allocTempArray<DTYPE>(&NAME##_device[_cudaGetCurrentDevice()],        \
                               STRINGIFY(NAME), {__VA_ARGS__}) != 0)           \
      return -1;                                                               \
    NAME = NAME##_device[_cudaGetCurrentDevice()];                             \
  }

#define PER_DEVICE_ARRAY_TEMP_DESTROY(NAME)                                    \
  {                                                                            \
    int device = _cudaGetCurrentDevice();                                      \
    _freeTempArray(NAME##_device[device]);                                     \
    NAME##_device[device] = nullptr;                                           \
  }

#define pow(x, y) powf(x, y)

/**
 * Miscellaneous helpers
 */

#define ENSURE_STREAM_PRIO(NAME, PRIORITY)                                     \
  static cudaStream_t NAME##_devs[MAX_CUDA_DEVICES] = {0};                     \
  /*[[maybe_unused]] cudaStream_t NAME = nullptr;*/                            \
  (void)NAME;                                                                  \
  {                                                                            \
    int device = _cudaGetCurrentDevice();                                      \
    if (device < 0)                                                            \
      return -1;                                                               \
    if (NAME##_devs[device] == nullptr) {                                      \
      int priority = PRIORITY;                                                 \
      CudaCheckedEx(cudaStreamCreateWithPriority(&NAME##_devs[device],         \
                                                 cudaStreamNonBlocking,        \
                                                 priority),                    \
                    "creating CUDA stream " STRINGIFY(NAME));                  \
    }                                                                          \
    NAME = NAME##_devs[device];                                                \
    cudaStreamSynchronize(NAME);                                               \
  };

#define ENSURE_STREAM(NAME) ENSURE_STREAM_PRIO(NAME, 0)

#define ENSURE_CACHE()                                                         \
  static jitify::JitCache caches[MAX_CUDA_DEVICES];                            \
  auto &kernel_cache = caches[0];

#define SET_DEVICE(dev)                                                        \
  {                                                                            \
    int cur = _cudaGetCurrentDevice();                                         \
    if (cur != dev && dev != -1)                                               \
      warn("!!! CUDA issue: expected device %d, current device %d", dev, cur); \
    cudaSetDevice(dev);                                                        \
  }

/**
 * @brief Raise a Python error.
 *
 * @param format The format string
 * @param args Format string arguments
 */
template <typename... Args>
void acquire_gil_and_raise_error(const std::string &format, Args... args) {
  std::string message = string_format(format, std::forward<Args>(args)...);
  critical(message);
#ifndef OPERATOR_STANDALONE
  PyGILState_STATE gstate;
  gstate = PyGILState_Ensure();
  PyErr_Format(PyExc_RuntimeError, message.c_str());
  PyGILState_Release(gstate);
#endif
}

/**
 * Check for pending CUDA errors, and raise a Python error
 * if one is found
 */
inline bool _cudaChecked(cudaError_t err, const char *file, int line,
                         const char *extra = nullptr) {
  if (err != cudaSuccess) {
    acquire_gil_and_raise_error("!!! CUDA Error in operator: %s:%d %s%s", file,
                                line, (extra != nullptr) ? extra : "",
                                cudaGetErrorString(err));

    // Attempt to clear the current CUDA error
    cudaDeviceSynchronize();
    cudaGetLastError();

    return false;
  }

  return true;
}

inline bool _ncclChecked(ncclResult_t err, const char *file, int line,
                         const char *extra = nullptr) {
  // if (err == ncclInProgress) {
  //     ncclResult_t state = err;
  //     do {
  //         ncclCommGetAsyncError(comm, &state);
  //     } while (state == ncclInProgress);
  //     err = state;
  // }

  if (err != ncclSuccess) {
    acquire_gil_and_raise_error("!!! NCCL Error in operator: %s:%d %s%s", file,
                                line, (extra != nullptr) ? extra : "",
                                ncclGetErrorString(err));

    // Attempt to clear the current CUDA error
    cudaDeviceSynchronize();
    cudaGetLastError();

    return false;
  }

  return true;
}

#define CudaCheckLaunch(grid, block)                                           \
  _cudaCheckKernelLaunch(grid, block, __FILE__, __LINE__)
inline void _cudaCheckKernelLaunch(dim3 grid, dim3 block, const char *file,
                                   int line) {
  cudaError_t err = cudaPeekAtLastError();
  if (err != cudaSuccess) {
    critical("!!! CUDA Error after kernel launch: %s:%d %s", file, line,
             cudaGetErrorString(err));
  }
}

template <typename T> int _freeTempArrayData(T *array) {
  CudaChecked(cudaFree(array->device_data));
  CudaChecked(cudaFreeHost(array->data));
  CudaChecked(cudaFreeHost(array->size));
  return 0;
}

template <typename T> int _freeTempArray(T *array) {
  if (_freeTempArrayData(array) != 0)
    return -1;

  CudaChecked(cudaFreeHost(array));
  return 0;
}

template <typename TData, typename T>
inline int _allocTempArray(T **array_ptr, const char *name,
                           std::initializer_list<int> dims) {
  int device = 0;
  CudaChecked(cudaGetDevice(&device));

  if (*array_ptr == nullptr) {
    CudaChecked(cudaMallocHost((void **)array_ptr, sizeof(T)));
    memset((void *)(*array_ptr), 0, sizeof(T));
  }

  T *array = *array_ptr;
  array->element_size = sizeof(TData);
  array->rank = dims.size();
  size_t nbytes = array->element_size;
  std::vector<int> dim_vec(dims);

  for (auto d : dim_vec) {
    nbytes *= (size_t)d;
  }

  if (array->nbytes != nbytes) {
    debug("%sallocating %llu bytes for %s for temporary data on device %d",
          array->nbytes > 0 ? "re" : "", nbytes, name, device);
    if (_freeTempArrayData(array) != 0)
      return -1;
    CudaChecked(cudaMallocHost((void **)(&array->data), nbytes));
    CudaChecked(cudaMalloc((void **)(&array->device_data), nbytes));
    CudaChecked(cudaMemset(array->device_data, 0, nbytes));
    CudaChecked(
        cudaMallocHost((void **)(&array->size), sizeof(size_t) * array->rank));
    for (int i = 0; i < array->rank; i++) {
      array->size[i] = dim_vec[i];
    }
    array->nbytes = nbytes;
  }

  return 0;
}
template <typename T> T round_down(T value, T factor) {}

template <typename... Args>
inline int _launchKernel(const char *file, int line, const char *kname,
                         jitify::KernelInstantiation &kernel,
                         std::tuple<dim3, dim3> tune, cudaStream_t stream,
                         dim3 mins, dim3 maxs, Args... args) {
  dim3 threads = std::get<0>(tune);
  dim3 tb = dim3(threads.x * threads.y * threads.z);
  dim3 threads_sub = std::get<1>(tune);
  dim3 block_sizes = dim3(threads.x * threads_sub.x, threads.y * threads_sub.y,
                          threads.z * threads_sub.z);

  // Convert our minimum values to block-aligned offsets
  dim3 offsets = dim3((mins.x / block_sizes.x) * block_sizes.x,
                      (mins.y / block_sizes.y) * block_sizes.y,
                      (mins.z / block_sizes.z) * block_sizes.z);

  // Resize the grid based on the block sizes and offsets
  dim3 grid(
      (int)(fmaxf(1.,
                  ceil((float)(maxs.x - offsets.x) / (float)(block_sizes.x)))),
      (int)(fmaxf(1.,
                  ceil((float)(maxs.y - offsets.y) / (float)(block_sizes.y)))),
      (int)(fmaxf(1,
                  ceil((float)(maxs.z - offsets.z) / (float)(block_sizes.z)))));

  if (grid.x >= 1 && grid.y >= 1 && grid.z >= 1) {
    // printf("%s: grid(%d, %d, %d), threads(%d, %d, %d), min=(%d, %d, %d),
    // max=(%d, %d, %d)\n", kname, grid.x, grid.y, grid.z, threads.x, threads.y,
    // threads.z, mins.x, mins.y, mins.z, maxs.x, maxs.y, maxs.z);
    _cudaChecked(
        (cudaError_t)(kernel.configure(grid, tb, 0, stream)
                          .launch(offsets, std::forward<Args>(args)...)),
        file, line);
  }

  return 0;
}

#define launchKernel(KNAME, ...)                                               \
  {                                                                            \
    if (_launchKernel(__FILE__, __LINE__, STRINGIFY(KNAME), KNAME##_tuned,     \
                      KNAME##_tune, __VA_ARGS__) < 0) {                        \
      return -1;                                                               \
    }                                                                          \
  }

inline uint64_t next_pow2(uint64_t x) {
  return (__builtin_popcount(x) == 1 || x == 1)
             ? x
             : 1 << (64 - __builtin_clzl(x - 1));
}

#define setupGrid(GRID, TB, THREAD, X, Y, Z)                                   \
  _setupGrid(#GRID, GRID, TB, THREAD, X, Y, Z)

#define setupGrid_v1(GRID, THREAD, X, Y, Z)                                    \
  {                                                                            \
    THREAD = dim3(1, 1, 1);                                                    \
    dim3 tmp;                                                                  \
    _setupGrid(#GRID, GRID, tmp, THREAD, X, Y, Z);                             \
  }
// Heuristic thread block sizing based on the grid
// Not perfect, but it'll do for now
inline void _setupGrid(const char *gridName, dim3 &grid, dim3 &cudaTb,
                       dim3 &threadBlock, int x_size, int y_size, int z_size) {
  if (threadBlock.x == threadBlock.y == threadBlock.z == 1) {
    if (z_size > 128) {
      threadBlock = dim3(1, 1, 64);
    } else if (z_size > 64) {
      threadBlock = dim3(1, 1, 32);
    } else if (y_size >= 8) {
      int y = std::max((int)next_pow2(std::min(128 / y_size, y_size)), 1);
      int z =
          std::max(std::min(64, (int)next_pow2(std::min(128 / y, z_size))), 1);
      threadBlock = dim3(1, y, z);
    } else if (x_size > 128) {
      threadBlock = dim3(128, 1, 1);
    } else {
      threadBlock = dim3(64, 1, 1);
    }
  }

  cudaTb = dim3(threadBlock.x * threadBlock.y * threadBlock.z, 1, 1);
  grid = dim3((int)ceil((float)x_size / (float)threadBlock.x),
              (int)ceil((float)y_size / (float)threadBlock.y),
              (int)ceil((float)z_size / (float)threadBlock.z));
}

inline int _cudaGetCurrentDevice() {
  int device = -1;
  CudaChecked(cudaGetDevice(&device));
  return device;
}

template <typename T> inline bool _cudaPointerIsAccessible(T *ptr) {
  struct cudaPointerAttributes attr = cudaPointerAttributes{};
  CudaChecked(cudaPointerGetAttributes(&attr, (const void *)ptr));
  return attr.devicePointer != NULL;
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
  return ret == cudaSuccess && attr.devicePointer != NULL;
}

#define _setCudaConstant(TYPE, NAME, VAL)                                      \
  {                                                                            \
    TYPE tmp = VAL;                                                            \
    CudaChecked(cudaMemcpyToSymbol(NAME, &tmp, sizeof(TYPE)));                 \
  }

template <typename T> size_t nbytes(const T *const __restrict obj) {
  size_t size = obj->element_size;
  for (int i = 0; i < obj->rank; i++) {
    size *= obj->size[i];
  }
  return size;
}

template <typename T>
int transferDataObject(cudaMemcpyKind kind, T *obj, size_t size = 0,
                       bool cond = true, cudaStream_t stream = nullptr,
                       const char *name = nullptr) {
  T *src = nullptr;
  T *dst = nullptr;
  if (kind == cudaMemcpyHostToDevice) {
    src = (T *)obj->data;
    dst = (T *)obj->device_data;
  } else {
    src = (T *)obj->device_data;
    dst = (T *)obj->data;
  }

  if (size == 0) {
    size = nbytes<T>(obj);
  }

  if (cond) {
    if (name != nullptr) {
      debug("transferring %s %s (%lx -> %lx)", name,
            kind == cudaMemcpyHostToDevice ? "H->D" : "D->H", src, dst);
    }
    if (stream != nullptr) {
      CudaChecked(cudaMemcpyAsync(dst, src, size, kind, stream));
    } else {
      CudaChecked(cudaMemcpy(dst, src, size, kind));
    }
  }

  return 0;
}

#define prepareDataObject(NAME, ...)                                           \
  _prepareDataObject(NAME, STRINGIFY(NAME), "prepareDataObject(" #NAME ")",    \
                     __VA_ARGS__)
template <typename T>
int _prepareDataObject(T *obj, const char *name, const char *nvtxRange,
                       size_t size = 0, bool copyIn = true,
                       cudaStream_t stream = nullptr) {
  nvtxRangePush(nvtxRange);
  int device = 0;
  int ret = 0;

  if (size == 0) {
    size = nbytes<T>(obj);
  }

  CudaChecked(cudaGetDevice(&device));
  if (!_cudaPtrIsManaged(obj->data)) {
    if (obj->device_data == nullptr) {
      debug("allocating %llu bytes for %s on device %d", size, name, device);
      CudaChecked(cudaMallocAsync((void **)&obj->device_data, size, stream));
      obj->operator_allocated = 1;
    } else {
      cudaPointerAttributes attrs = cudaPointerAttributes{};
      CudaChecked(cudaPointerGetAttributes(&attrs, obj->device_data));

      if (attrs.devicePointer == NULL) {
        acquire_gil_and_raise_error(
            "!!! CUDA error: object %s with host ptr %llx,"
            "is inaccessible by device %d (created on device %d)",
            name, obj->data, device, attrs.device);
        nvtxRangePop();
        return -1;
      } else if (attrs.device != device) {
        critical(
            "!!! CUDA issue: object %s with host ptr %llx, original device "
            "pointer %llx, local device pointer %llx was created on device %d "
            "and is being accessed by device %d - this should not happen",
            name, obj->data, obj->device_data, attrs.devicePointer,
            attrs.device, device);
      }
    }
    ret = transferDataObject(cudaMemcpyHostToDevice, obj, size, copyIn, stream,
                             name);
  }
  nvtxRangePop();
  return ret;
}
#define destroyDataObject(NAME, ...)                                           \
  {                                                                            \
    if (_destroyDataObject(NAME, STRINGIFY(NAME), __VA_ARGS__) < 0)            \
      return -1;                                                               \
  }

template <typename T>
int _destroyDataObject(T *obj, const char *name, bool del = true,
                       cudaStream_t stream = nullptr) {
  if (del && obj->operator_allocated) {
    debug("freeing %s", name);
    CudaChecked(cudaFreeAsync(obj->device_data, stream));
    obj->device_data = nullptr;
  }

  return 0;
}

using tuned_kernel = std::pair<dim3, dim3>;
typedef std::map<std::shared_ptr<jitify::detail::CUDAKernel>, tuned_kernel>
    tuningDict;

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

static float _occupancyForKernel(CUfunction &k, const dim3 &block) {
  int max_sm_registers = 0;
  int max_block_registers;
  int max_sm_threads;
  int max_sm_blocks;

  int device = 0;
  cudaGetDevice(&device);
  cudaDeviceGetAttribute(&max_sm_threads,
                         cudaDevAttrMaxThreadsPerMultiProcessor, device);
  cudaDeviceGetAttribute(&max_sm_blocks, cudaDevAttrMaxBlocksPerMultiprocessor,
                         device);
  cudaDeviceGetAttribute(&max_sm_registers,
                         cudaDevAttrMaxRegistersPerMultiprocessor, device);
  cudaDeviceGetAttribute(&max_block_registers,
                         cudaDevAttrMaxRegistersPerMultiprocessor, device);

  int regs;
  cuFuncGetAttribute(&regs, CU_FUNC_ATTRIBUTE_NUM_REGS, k);

  int block_regs = regs * block.x * block.y * block.z;
  int warp_regs = regs * 32;

  int warps = (int)ceil((block.x * block.y * block.z) / 32.f);

  int max_sm_warps = max_sm_threads / 32;

  int max_warps_per_sm_reg = (int)ceil(max_sm_registers / (float)warp_regs);
  int max_block_per_sm_reg = (int)ceil(max_sm_registers / (float)block_regs);

  int active_warps =
      std::min(max_warps_per_sm_reg, max_block_per_sm_reg * warps);

  int max_block_per_sm_warp = (int)ceil(max_sm_warps / (float)warps);

  debug("warps per sm (register limited) = %d", active_warps);
  debug("max blocks per sm (register limited) = %d", max_block_per_sm_reg);
  debug("max blocks per sm (thread limited) = %d", max_block_per_sm_warp);

  float warp_sm_reg_occupancy =
      (32.0f * (float)active_warps) / (float)max_sm_threads;
  float warp_sm_occupancy =
      (float)(max_block_per_sm_warp * block.x * block.y * block.z) /
      (float)max_sm_threads;
  float block_sm_reg_occupancy =
      ((block.x * block.y * block.z) * (float)max_block_per_sm_reg) /
      (float)max_sm_blocks;
  return fminf(1.0f, fminf(warp_sm_reg_occupancy,
                           fminf(warp_sm_occupancy, block_sm_reg_occupancy)));
}

std::mutex kernel_compile_mutex;

static bool
_check_kernel(const dim3 &block, const dim3 &sub_block,
              std::function<jitify::KernelInstantiation(dim3, dim3)> &builder,
              bool &is_valid, float &est_occupancy, int &max_block, int &regs,
              float &occupancy) {
  // verify that it works by checking occupancy with the new block size
  // Keep kernel_inst alive to ensure the CUfunction remains valid
  jitify::KernelInstantiation kernel_inst;
  CUfunction k = NULL;
  {
    std::lock_guard<std::mutex> lock(kernel_compile_mutex);
    kernel_inst = builder(block, sub_block);
    k = kernel_inst.get_function();
    cuFuncSetCacheConfig(k, CU_FUNC_CACHE_PREFER_L1);
  }
  int grid = 0;

  is_valid = false;
  est_occupancy = 0.f;
  regs = 0;

  debug("trying with block size (%d, %d, %d)..", block.x, block.y, block.z);

  CUresult res =
      cuOccupancyMaxPotentialBlockSize(&grid, &max_block, k, nullptr, 0, 0);
  if (res == 0) {

    if (cuFuncGetAttribute(&regs, CU_FUNC_ATTRIBUTE_NUM_REGS, k) != 0)
      return false;

    occupancy = _occupancyForKernel(k, block);
    debug("rebuilt kernel is valid, has max occupancy at %d threads "
          "(%.2f device occupancy), uses %d registers",
          max_block, 100. * occupancy, regs);
    if (max_block >= block.x * block.y * block.z) {
      is_valid = true;

      est_occupancy = (float)(block.x * block.y * block.z) / max_block;

      return true;
    }
  }

  return false;
}

static float _est_waste(dim3 grid, dim3 block) {
  long block_pt = block.x * block.y * block.z;
  long grid_pt = grid.x * grid.y * grid.z;
  return (float)(grid_pt % block_pt) / (float)grid_pt;
}

inline float _est_efficiency(int ideal, dim3 proposed_block) {
  return (float)(proposed_block.x * proposed_block.y * proposed_block.z) /
         (float)ideal;
}

inline bool compare_options(float occupancy1, const dim3 &block1,
                            float occupancy2, const dim3 &block2) {
  float diff = fabs(occupancy2 - occupancy1);
  if (occupancy1 > occupancy2 && diff > 0.01)
    return false;

  if (diff < 0.01)
    return (block1.x * block1.y * block1.z) < (block2.x * block2.y * block2.z);

  return true;
}

static tuned_kernel
performTuning(tuningDict &tuning, const char *name, dim3 preferred,
              dim3 preferred_sub, dim3 expected_grid, int max_block_dimension,
              std::function<jitify::KernelInstantiation(dim3, dim3)> builder) {
  nvtxRangePush(name);
  int device = 0;
  int max_sm_resident_blocks = 0;
  int sm_count = 0;

  // Keep kernel_inst alive to ensure CUfunction remains valid throughout tuning
  jitify::KernelInstantiation kernel_inst;
  std::shared_ptr<jitify::detail::CUDAKernel> kernel_ptr;
  CUfunction cf = NULL;

  cudaGetDevice(&device);
  cudaDeviceGetAttribute(&max_sm_resident_blocks,
                         cudaDevAttrMaxBlocksPerMultiprocessor, device);
  cudaDeviceGetAttribute(&sm_count, cudaDevAttrMultiProcessorCount, device);
  {
    std::lock_guard<std::mutex> lock(kernel_compile_mutex);

    kernel_inst = builder(preferred, preferred_sub);
    kernel_ptr = kernel_inst.cuda_kernel_ptr();

    // Calculate the maximum possible occupancy for the preferred block size
    cf = kernel_inst.get_function();
    if (tuning.find(kernel_ptr) != tuning.end()) {
      nvtxRangePop();
      return tuning.at(kernel_ptr);
    }
  }

  for (int i = max_block_dimension; i < 3; i++) {
    dim3_set(preferred, i, 1);
  }

  debug("=== performing tuning for %s with default block "
        "%d,%d,%d, sub-block "
        "%d,%d,%d, and expected grid size (%d, %d, %d)",
        name, preferred.x, preferred.y, preferred.z, preferred_sub.x,
        preferred_sub.y, preferred_sub.z, expected_grid.x, expected_grid.y,
        expected_grid.z);

  int grid;
  int max_block = 0;

  CUresult res =
      cuOccupancyMaxPotentialBlockSize(&grid, &max_block, cf, nullptr, 0, 0);

  if (res != 0) {
    critical("!!! invalid kernel detected! could not calculate occupancy "
             "for %s "
             "with blocksize (%d, %d, %d), sub-block size (%d, %d, %d)!",
             name, preferred.x, preferred.y, preferred.z, preferred_sub.x,
             preferred_sub.y, preferred_sub.z);
    nvtxRangePop();
    assert(false);
  }

  debug("maximum occupancy is at %d threads with minimum %d grid size",
        max_block, grid);

  int threads = max_block;

  int last_dim = dim3_get(expected_grid, max_block_dimension - 1);

  int warps = (int)ceil(threads / 32.f);

  debug("maximum warps = %d", warps);
  assert(warps >= 1);

  tuned_kernel result(preferred, preferred_sub);

  dim3 block(32, 1, 1);

  // Check for small dimensions
  bool small_dims[3] = {false};
  bool has_small_dims = false;
  for (int d = 0; d < max_block_dimension; d++) {
    int dv = dim3_get(expected_grid, d);
    if (dv < 32 || ((dv % 32) / dv) > 0.75) {
      small_dims[d] = true;
      has_small_dims = true;
    }
  }

  if (has_small_dims)
    debug("expected grid (%d, %d, %d) has one or more small dimensions "
          "(<32 or >75%% expected waste)",
          expected_grid.x, expected_grid.y, expected_grid.z);

  int nearest_square = 0;

  // Assume the initial input is 'best' until we know otherwise
  float best_eff = _occupancyForKernel(cf, preferred);
  debug("base occupancy is %.2f", best_eff);
  if (best_eff > 0.66) {
    std::lock_guard<std::mutex> lock(kernel_compile_mutex);
    tuning[kernel_ptr] = result;
    nvtxRangePop();
    return tuning[kernel_ptr];
  }

  float next_eff = 0.f;
  bool next_valid = false;
  float tmp = 0.f;

  dim3 next_block;

  std::set<std::tuple<int, int, int>> tried;
  int next_regs;

  switch (max_block_dimension) {
  case 1:
    // We're done, 1D is a simple case
    debug("kernel is 1D and maximum block size is %d", max_block);
    std::get<0>(result).x = max_block;
    break;

  case 2:
    // for now, just use preferred size, even though that's suboptimal
    // nearly always

    break;

  case 3:
    // Short-circuit if we're already using less blocks than the device can
    // simultaneously run
    if (best_eff > 0.1 &&
        sm_count * max_sm_resident_blocks >
            (expected_grid.x * expected_grid.y * expected_grid.z) /
                (block.x * block.y * block.z)) {
      debug("Required block count ~%d already fits on device (%d SM, %d "
            "block/SM). Giving up now.",
            (expected_grid.x * expected_grid.y * expected_grid.z) /
                (block.x * block.y * block.z),
            sm_count, max_sm_resident_blocks);
      break;
    }
    // Force small dimensions to have small block sizes, even if it means
    // warp divergence
    for (int d = 0; d < 3; d++) {
      if (small_dims[d])
        dim3_set(block, d, 4);
    }

    // try and make a better guess
    nearest_square =
        std::max(1, (int)(floor(sqrtf(std::max(1, threads) / (float)block.x))));
    if ((small_dims[0] || small_dims[1]) &&
        _est_waste(expected_grid,
                   dim3(nearest_square, nearest_square, block.z)) > 0.2) {
      debug("remaining small dimensions would lead to excess wasted "
            "grid points, not trying a square");
    } else {
      next_block = dim3(block.x,
                        block.y > 1 ? std::min((int)block.y, nearest_square)
                                    : nearest_square,
                        block.z > 1 ? std::min((int)block.z, nearest_square)
                                    : nearest_square);
      _check_kernel(next_block, preferred_sub, builder, next_valid, tmp,
                    max_block, next_regs, next_eff);

      if (next_valid) {
        debug("nearest square attempt would be %d,%d,%d "
              "(occupancy %.2f%%)",
              next_block.x, next_block.y, next_block.z, 100. * next_eff);

        if (compare_options(best_eff, std::get<0>(result), next_eff,
                            next_block)) {
          result = tuned_kernel(next_block, preferred_sub);
          best_eff = next_eff;
        }
      }

      tried.emplace(std::make_tuple(next_block.x, next_block.y, next_block.z));

      int yz = std::max(1, (int)ceil((float)max_block / (float)block.x));

      int x = (int)(floor(sqrtf(yz)));
      while (x > 1 && (yz % x > 0))
        x--;

      // find largest grid for which (kernel max occupancy threads) %
      // (x*y*z) ~= 0
      for (int offset = 0; offset < 3 && best_eff < 0.99999; offset++) {
        int next_y = x - offset;
        if (next_y <= 0 || yz % next_y != 0)
          continue;
        next_block = dim3(block.x, (yz / next_y), next_y);

        if (tried.find(std::make_tuple(next_block.x, next_block.y,
                                       next_block.z)) == tried.end()) {
          debug("next guess is %d,%d,%d", next_block.x, next_block.y,
                next_block.z);
          _check_kernel(next_block, preferred_sub, builder, next_valid, tmp,
                        max_block, next_regs, next_eff);

          debug(" (occupancy %.2f%%)", 100. * next_eff);
          if (next_valid && compare_options(best_eff, std::get<0>(result),
                                            next_eff, next_block)) {
            best_eff = next_eff;
            result = tuned_kernel(next_block, preferred_sub);
          }
          tried.emplace(
              std::make_tuple(next_block.x, next_block.y, next_block.z));
        }
      }
    }

    break;
  default:
    break;
  }

  debug("selected block (%d, %d, %d), sub_block (%d, %d, %d) with est. "
        "occupancy = %.2f",
        std::get<0>(result).x, std::get<0>(result).y, std::get<0>(result).z,
        std::get<1>(result).x, std::get<1>(result).y, std::get<1>(result).z,
        best_eff);
  {
    std::lock_guard<std::mutex> lock(kernel_compile_mutex);
    tuning[kernel_ptr] = result;
  }
  nvtxRangePop();
  return result;
}

static bool exceptionOccured() {
#ifndef OPERATOR_STANDALONE
  PyGILState_STATE gstate;
  gstate = PyGILState_Ensure();
  bool exceptionOccured = PyErr_Occurred() != NULL;
  PyGILState_Release(gstate);
  return exceptionOccured;
#else
  return false;
#endif
}

template <typename TDataobj>
static int devito_cuda_async_gather_4d(float *__restrict buf, TDataobj dataobj,
                                       const int x_sz, const int y_sz,
                                       const int z_sz, const int w_ofs,
                                       const int x_ofs, const int y_ofs,
                                       const int z_ofs, cudaStream_t stream) {
  struct cudaMemcpy3DParms copy_params = {0};
  assert(dataobj->rank == 4);

  int rank = dataobj->rank;
  size_t pitch = dataobj->size[rank - 1] * sizeof(float);
  size_t row_size = dataobj->size[rank - 1] * sizeof(float);
  size_t ptr_offset = w_ofs * dataobj->size[rank - 3] *
                      dataobj->size[rank - 2] * dataobj->size[rank - 1];

  copy_params.srcPtr =
      make_cudaPitchedPtr(((float *)(dataobj->device_data)) + ptr_offset, pitch,
                          row_size, dataobj->size[rank - 2]);
  copy_params.srcPos = make_cudaPos(sizeof(float) * (z_ofs), (y_ofs), (x_ofs));
  copy_params.dstPtr = make_cudaPitchedPtr(buf, sizeof(float) * (z_sz),
                                           sizeof(float) * (z_sz), y_sz);
  copy_params.dstPos = make_cudaPos(sizeof(float) * (0), (0), (0));
  copy_params.extent = make_cudaExtent(sizeof(float) * (z_sz), y_sz, x_sz);
  copy_params.kind = cudaMemcpyDeviceToDevice;
  CudaChecked(cudaMemcpy3DAsync(&copy_params, stream));
  return cudaSuccess;
}

template <typename TDataobj>
static int devito_cuda_async_scatter_4d(TDataobj dataobj,
                                        const float *const __restrict buf,
                                        const int x_sz, const int y_sz,
                                        const int z_sz, const int w_ofs,
                                        const int x_ofs, const int y_ofs,
                                        const int z_ofs, cudaStream_t stream) {
  struct cudaMemcpy3DParms copy_params = {0};
  assert(dataobj->rank == 4);
  int rank = dataobj->rank;
  size_t pitch = dataobj->size[rank - 1] * sizeof(float);
  size_t row_size = dataobj->size[rank - 1] * sizeof(float);
  size_t ptr_offset = w_ofs * dataobj->size[rank - 3] *
                      dataobj->size[rank - 2] * dataobj->size[rank - 1];
  copy_params.dstPtr =
      make_cudaPitchedPtr(((float *)(dataobj->device_data)) + ptr_offset, pitch,
                          row_size, dataobj->size[rank - 2]);
  copy_params.dstPos = make_cudaPos(sizeof(float) * (z_ofs), (y_ofs), (x_ofs));
  copy_params.srcPtr = make_cudaPitchedPtr((float *)buf, sizeof(float) * (z_sz),
                                           sizeof(float) * (z_sz), y_sz);
  copy_params.srcPos = make_cudaPos(sizeof(float) * (0), (0), (0));
  copy_params.extent = make_cudaExtent(sizeof(float) * (z_sz), y_sz, x_sz);
  copy_params.kind = cudaMemcpyDeviceToDevice;
  CudaChecked(cudaMemcpy3DAsync(&copy_params, stream));
  return cudaSuccess;
}

// Multi-function halo update using NCCL
//
// Basically just exists to avoid Devito spitting out a ton of nearly
// identical halo exchange functions.
//
// That, and we can do C++-y things more easily in a C++ header than in
// Devito's AST/IR format.
//
template <typename TDataobj, typename TMPIMsg>
static int devito_cuda_async_multi_haloupdate(
    std::initializer_list<TDataobj> functions,
    std::initializer_list<TMPIMsg> msgs, int otime, int ncomms,
    cudaStream_t nccl_stream, cudaStream_t kernel_stream,
    ncclComm_t nccl_comm) {
  assert(functions.size() == msgs.size());
  cudaEvent_t update_ev = nullptr;
  CudaChecked(cudaEventCreateWithFlags(&update_ev, cudaEventDisableTiming));
  for (ptrdiff_t i = 0; i < ncomms; i++) {
    for (ptrdiff_t f = 0; f < functions.size(); f++) {
      auto &msg = *(msgs.begin() + f);
      auto &function = *(functions.begin() + f);
      if (msg[i].torank != MPI_PROC_NULL) {
        Checked(devito_cuda_async_gather_4d<TDataobj>(
            (float *)msg[i].bufg, function, msg[i].sizes[0], msg[i].sizes[1],
            msg[i].sizes[2], otime, msg[i].ofsg[0], msg[i].ofsg[1],
            msg[i].ofsg[2], kernel_stream));
      }
    }
  }

  CudaChecked(cudaEventRecord(update_ev, kernel_stream));
  CudaChecked(cudaStreamWaitEvent(nccl_stream, update_ev, 0));
  CudaChecked(cudaEventDestroy(update_ev));
  NcclChecked(ncclGroupStart());
  for (ptrdiff_t i = 0; i < ncomms; i++) {
    for (ptrdiff_t f = 0; f < functions.size(); f++) {
      auto &msg = *(msgs.begin() + f);
      if (msg[i].fromrank != MPI_PROC_NULL) {
        NcclChecked(ncclRecv(
            msg[i].bufs, msg[i].sizes[0] * msg[i].sizes[1] * msg[i].sizes[2],
            ncclFloat32, msg[i].fromrank, nccl_comm, nccl_stream));
      }

      if (msg[i].torank != MPI_PROC_NULL) {
        NcclChecked(ncclSend(
            msg[i].bufg, msg[i].sizes[0] * msg[i].sizes[1] * msg[i].sizes[2],
            ncclFloat32, msg[i].torank, nccl_comm, nccl_stream));
      }
    }
  }
  NcclChecked(ncclGroupEnd());
  return 0;
}

template <typename TDataobj, typename TMPIMsg>
static int
devito_cuda_async_multi_halowait(std::initializer_list<TDataobj> functions,
                                 std::initializer_list<TMPIMsg> msgs, int otime,
                                 int ncomms, cudaStream_t nccl_stream,
                                 cudaStream_t kernel_stream,
                                 ncclComm_t nccl_comm) {
  assert(functions.size() == msgs.size());
  cudaEvent_t update_ev = nullptr;
  CudaChecked(cudaEventCreateWithFlags(&update_ev, cudaEventDisableTiming));
  CudaChecked(cudaEventRecord(update_ev, nccl_stream));
  CudaChecked(cudaStreamWaitEvent(kernel_stream, update_ev));

  for (ptrdiff_t i = 0; i < ncomms; i++) {
    for (ptrdiff_t f = 0; f < functions.size(); f++) {
      auto &msg = *(msgs.begin() + f);
      auto &function = *(functions.begin() + f);
      if (msg[i].fromrank != MPI_PROC_NULL) {
        Checked(devito_cuda_async_scatter_4d<TDataobj>(
            function, (float *)msg[i].bufg, msg[i].sizes[0], msg[i].sizes[1],
            msg[i].sizes[2], otime, msg[i].ofsg[0], msg[i].ofsg[1],
            msg[i].ofsg[2], kernel_stream));
      }
    }
  }

  return 0;
}

namespace pair_iterators {
template <typename T1, typename T2>
std::pair<T1, T2> operator++(std::pair<T1, T2> &it) {
  ++it.first;
  ++it.second;
  return it;
}
} // namespace pair_iterators

template <typename TDataobj>
static int
devito_cuda_async_d2h_destroy_many(std::initializer_list<TDataobj> dataobjs,
                                   std::initializer_list<char *> names,
                                   bool devicerm, bool updatehost,
                                   cudaStream_t stream) {
  assert(dataobjs.size() == names.size());
  for (auto its = std::make_pair(dataobjs.begin(), names.begin()),
            end = std::make_pair(dataobjs.end(), names.end());
       its != end; ++its) {
    if (transferDataObject(cudaMemcpyDeviceToHost, *(its.first), 0, updatehost,
                           stream, *(its.second)) < 0 ||
        _destroyDataObject(*(its.first), *(its.second), devicerm, stream) < 0) {
      return -1;
    }
  }
  return 0;
}

template <typename TDataobj>
static int
devito_cuda_async_h2d_prepare_many(std::initializer_list<TDataobj> dataobjs,
                                   std::initializer_list<char *> names,
                                   bool devicecreate, bool updatedevice,
                                   cudaStream_t stream) {
  assert(dataobjs.size() == names.size());
  char nvtxRange[256];
  for (auto its = std::make_pair(dataobjs.begin(), names.begin()),
            end = std::make_pair(dataobjs.end(), names.end());
       its != end; ++its) {
    memset(nvtxRange, 0, 256);
    snprintf(nvtxRange, 256, "prepareDataObject(%s)", *(its.second));
    if (_prepareDataObject(*(its.first), *(its.second), nvtxRange, 0,
                           devicecreate || updatedevice, stream) < 0) {
      return -1;
    }
  }
  return 0;
}

class CudaSectionTimer {
public:
  CudaSectionTimer(cudaStream_t stream, double *section_ptr)
      : _stream(stream), _section_ptr(section_ptr), _start_event(nullptr),
        _end_event(nullptr) {}

  void start() {
    if (_start_event) {
      return;
    }
    cudaEventCreate(&_start_event);
    cudaEventRecord(_start_event, _stream);
  }

  void stop() {
    cudaEventCreate(&_end_event);
    cudaEventRecord(_end_event, _stream);
  }

  ~CudaSectionTimer() {
    if (_start_event) {
      cudaEventDestroy(_start_event);
      _start_event = nullptr;
    }
    if (_end_event) {
      cudaEventDestroy(_end_event);
      _end_event = nullptr;
    }
  }

  void resolve() {
    if (_start_event == nullptr) {
      return;
    }

    if (_end_event == nullptr) {
      debug("CudaSectionTimer: stop() was not called, possible codegen bug");
      return;
    }
    cudaEventSynchronize(_end_event);

    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, _start_event, _end_event);

    *_section_ptr +=
        static_cast<double>(milliseconds) / 1000.0; // Convert to seconds

    cudaEventDestroy(_start_event);
    cudaEventDestroy(_end_event);
    _start_event = nullptr;
    _end_event = nullptr;
  }

private:
  cudaStream_t _stream;
  double *_section_ptr;
  cudaEvent_t _start_event;
  cudaEvent_t _end_event;
};

class CudaSectionTimers {
public:
  CudaSectionTimers(cudaStream_t timer_stream) : _timer_stream(timer_stream) {
    // This should probably be enough - normally we target a few thousand time
    // steps at most, and large elastic operators may have 30-50 sections.
    _timers.reserve(131072);
  }

  ~CudaSectionTimers() { _timers.clear(); }

  CudaSectionTimer &startNewTimer(double *section_ptr) {
    _timers.emplace_back(_timer_stream, section_ptr);
    _timers.back().start();
    return _timers.back();
  }

  void resolveTimers() {
    for (auto &timer : _timers) {
      timer.resolve();
    }
  }

private:
  std::vector<CudaSectionTimer> _timers;
  cudaStream_t _timer_stream;
};

#define DEVITO_CUDA_PROLOGUE()                                                 \
  CudaSectionTimers _cuda_section_timers(kernel_stream);                       \
  nvtxRangePush(__FUNCTION__);

#define CUDA_START_TIMER(T, S)                                                 \
  auto &_timer_##S = _cuda_section_timers.startNewTimer(&(T->S));

#define CUDA_STOP_TIMER(ST) _timer_##ST.stop();

#define DEVITO_CUDA_EPILOGUE()                                                 \
  do {                                                                         \
    _cuda_section_timers.resolveTimers();                                      \
    if (devicerm || updatehost) {                                              \
      CudaChecked(cudaStreamSynchronize(kernel_stream));                       \
      nvtxRangePop();                                                          \
    }                                                                          \
  } while (0);

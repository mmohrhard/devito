#pragma once

#include "nvtx3/nvToolsExt.h"
#include <cmath>
#include <cuda.h>
#include <devito/jitify.hpp>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <tuple>
#include <utility>

#define STRINGIFY(x) _stringify(x)
#define _stringify(x) #x

#ifndef NVRTC_CUDA_ARCH
// Default to compute capability 7.0 aka V100
#define NVRTC_CUDA_ARCH compute_70
#endif

#ifdef KERNEL_DEBUGGING
#define DEBUG_OPTS "-G",
#else
#define DEBUG_OPTS
#endif

#define NVRTC_OPTS                                                             \
  {                                                                            \
    "--ftz=true", "--fmad=true", "--prec-sqrt=true", "--prec-div=true",        \
        DEBUG_OPTS "--modify-stack-limit=false", "--split-compile=4",          \
        "--gpu-architecture=" STRINGIFY(NVIDIA_CUDA_ARCH),                     \
        "--extra-device-vectorization", "--minimal", "--restrict"                            \
  }

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
  if (_logHandler == nullptr)
    return;

  std::string message = string_format(format, std::forward<Args>(args)...);
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

#define CudaChecked(f) _cudaChecked((f), __FILE__, __LINE__)
#define CudaCheckedEx(f, msg) _cudaChecked((f), __FILE__, __LINE__, msg)

#define MAX_CUDA_DEVICES 16
#define PER_DEVICE_TEMP_GET(T, NAME, SIZE)                                     \
  T *NAME = nullptr;                                                           \
  static T *NAME##_device[MAX_CUDA_DEVICES] = {0};                             \
  {                                                                            \
    int device = 0;                                                            \
    cudaGetDevice(&device);                                                    \
    assert(device >= 0 && device < MAX_CUDA_DEVICES);                          \
    if (NAME##_device[device] == nullptr) {                                    \
      debug("allocating %llu bytes for " STRINGIFY(NAME) " on device %d",      \
            SIZE, device);                                                     \
      CudaChecked(cudaMallocAsync((void **)&NAME##_device[device], (SIZE),     \
                                  kernel_stream));                             \
      cudaMemset(NAME##_device[device], 1, (SIZE));                            \
    }                                                                          \
    NAME = NAME##_device[device];                                              \
  }

#define PER_DEVICE_TEMP_DESTROY(NAME)                                          \
  {                                                                            \
    int device = 0;                                                            \
    cudaGetDevice(&device);                                                    \
    assert(device >= 0 && device < MAX_CUDA_DEVICES);                          \
    CudaChecked(cudaFreeAsync(NAME##_device[device], kernel_stream));          \
    NAME##_device[device] = nullptr;                                           \
  }

#define PER_DEVICE_ARRAY_TEMP_DECLARE(NAME, ARRAYTYPE)                         \
  static ARRAYTYPE *NAME##_device[MAX_CUDA_DEVICES] = {0};

#define PER_DEVICE_ARRAY_TEMP_GET(NAME, NBYTES)                                \
  _allocTempArray(&NAME##_device[_cudaGetCurrentDevice()], NBYTES,             \
                  STRINGIFY(NAME))

#define PER_DEVICE_ARRAY_TEMP_DESTROY(NAME)                                    \
  {                                                                            \
    int device = _cudaGetCurrentDevice();                                      \
    _freeTempArray(NAME##_device[device]);                                     \
    NAME##_device[device] = nullptr;                                           \
  }

#define pow(x, y) powf(x, y)

#define ENSURE_STREAM(NAME)                                                    \
  static cudaStream_t NAME##_devs[MAX_CUDA_DEVICES] = {0};                     \
  [[maybe_unused]] cudaStream_t NAME = nullptr;                                \
  {                                                                            \
    int device = _cudaGetCurrentDevice();                                      \
    if (NAME##_devs[device] == nullptr)                                        \
      CudaCheckedEx(cudaStreamCreateWithFlags(&NAME##_devs[device],            \
                                              cudaStreamNonBlocking),          \
                    "creating CUDA stream " STRINGIFY(NAME));                  \
                                                                               \
    NAME = NAME##_devs[device];                                                \
  };

#define ENSURE_CACHE()                                                         \
  static jitify::JitCache caches[MAX_CUDA_DEVICES];                            \
  auto &kernel_cache = caches[0];

inline void _cudaChecked(cudaError_t err, const char *file, int line,
                         const char *extra = nullptr) {
  if (err != cudaSuccess) {
    err = cudaGetLastError();
    critical("!!! CUDA Error in operator: %s:%d %s%s", file, line,
             (extra != nullptr) ? extra : "", cudaGetErrorString(err));
  }
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

template <typename T> void _freeTempArrayData(T *array) {
  CudaChecked(cudaFree(array->device_data));
  CudaChecked(cudaFreeHost(array->data));
}
template <typename T> void _freeTempArray(T *array) {
  _freeTempArrayData(array);
  CudaChecked(cudaFreeHost(array));
}

template <typename T>
inline T *_allocTempArray(T **array_ptr, size_t nbytes, const char *name) {
  int device = 0;
  cudaGetDevice(&device);

  if (*array_ptr == nullptr) {
    CudaChecked(cudaMallocHost((void **)array_ptr, sizeof(T)));
    memset((void *)(*array_ptr), 0, sizeof(T));
  }

  T *array = *array_ptr;
  if (array->nbytes != nbytes) {
    debug("%sallocating %llu bytes for %s for temporary data on device %d",
          array->nbytes > 0 ? "re" : "", nbytes, name, device);
    _freeTempArrayData(array);
    CudaChecked(cudaMallocHost((void **)(&array->data), nbytes));
    CudaChecked(cudaMalloc((void **)(&array->device_data), nbytes));
    cudaMemset(array->device_data, 0, nbytes);
    array->nbytes = nbytes;
  }

  return array;
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
      int y = max((int)next_pow2(min(128 / y_size, y_size)), 1);
      int z = max(min(64, (int)next_pow2(min(128 / y, z_size))), 1);
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

template <typename T> inline void _cudaEnsureAllocated() {}

#define _setCudaConstant(TYPE, NAME, VAL)                                      \
  {                                                                            \
    TYPE tmp = VAL;                                                            \
    CudaChecked(cudaMemcpyToSymbol(NAME, &tmp, sizeof(TYPE)));                 \
  }

template <typename T>
void transferDataObject(cudaMemcpyKind kind, T *obj, size_t size = 0,
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

  if (size > 0 && cond) {
    if (stream != nullptr)
      CudaChecked(cudaMemcpyAsync(dst, src, size, kind, stream));
    else {
      debug("transferring %s %s (%lx -> %lx)", name,
            kind == cudaMemcpyHostToDevice ? "H->D" : "D->H", src, dst);

      CudaChecked(cudaMemcpy(dst, src, size, kind));
    }
  }
}

#define prepareDataObject(NAME, ...)                                           \
  _prepareDataObject(NAME, STRINGIFY(NAME), "prepareDataObject(" #NAME ")",    \
                     __VA_ARGS__);
template <typename T>
void _prepareDataObject(T *obj, const char *name, const char *nvtxRange,
                        size_t size, bool copyIn = true,
                        cudaStream_t stream = nullptr) {
  nvtxRangePush(nvtxRange);
  int device = 0;
  cudaGetDevice(&device);
  if (!_cudaPtrIsManaged(obj->data)) {
    if (obj->device_data == nullptr) {
      debug("allocating %llu bytes for %s on device %d", size, name, device);
      CudaChecked(cudaMallocAsync((void **)&obj->device_data, size, stream));
      obj->operator_allocated = 1;
    }
    transferDataObject(cudaMemcpyHostToDevice, obj, size, copyIn, stream, name);
  }
  nvtxRangePop();
}
#define destroyDataObject(NAME, ...)                                           \
  _destroyDataObject(NAME, STRINGIFY(NAME), __VA_ARGS__);

template <typename T>
void _destroyDataObject(T *obj, const char *name, bool del = true,
                        cudaStream_t stream = nullptr) {
  if (del && obj->operator_allocated) {
    debug("freeing %s", name);
    CudaChecked(cudaFreeAsync(obj->device_data, stream));
    obj->device_data = nullptr;
  }
}

using tuned_kernel = std::pair<dim3, dim3>;
typedef std::map<CUfunction, tuned_kernel> tuningDict;

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

  int warps = (block.x * block.y * block.z) / 32;

  int max_sm_warps = max_sm_threads / 32;

  int max_warps_per_sm_reg = max_sm_registers / warp_regs;
  int max_block_per_sm_reg = max_sm_registers / block_regs;

  int active_warps = min(max_warps_per_sm_reg, max_block_per_sm_reg * warps);

  int max_block_per_sm_warp = max_sm_warps / warps;

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

static bool
_check_kernel(const dim3 &block, const dim3 &sub_block,
              std::function<jitify::KernelInstantiation(dim3, dim3)> &builder,
              bool &is_valid, float &est_occupancy, int &max_block, int &regs,
              float &occupancy) {
  // verify that it works by checking occupancy with the new block size
  CUfunction k = builder(block, sub_block);

  int grid = 0;

  is_valid = false;
  est_occupancy = 0.f;
  regs = 0;

  debug("trying with block size (%d, %d, %d)..", block.x, block.y, block.z);

  CUresult res = cuOccupancyMaxPotentialBlockSize(&grid, &max_block,
                                                  (CUfunction)k, nullptr, 0, 0);
  if (res == 0) {

    if (cuFuncGetAttribute(&regs, CU_FUNC_ATTRIBUTE_NUM_REGS, (CUfunction)k) !=
        0)
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

std::mutex cache_mutex;

static tuned_kernel
performTuning(tuningDict &tuning, const char *name, dim3 preferred,
              dim3 preferred_sub, dim3 expected_grid, int max_block_dimension,
              std::function<jitify::KernelInstantiation(dim3, dim3)> builder) {
  nvtxRangePush(name);
  int device = 0;
  int max_sm_resident_blocks = 0;
  int sm_count = 0;

  cudaGetDevice(&device);
  cudaDeviceGetAttribute(&max_sm_resident_blocks,
                         cudaDevAttrMaxBlocksPerMultiprocessor, device);
  cudaDeviceGetAttribute(&sm_count, cudaDevAttrMultiProcessorCount, device);

  int max_block = 0;

  auto tmp_kernel = builder(preferred, preferred_sub);

  // Calculate the maximum possible occupancy for the preferred block size
  CUfunction cf = (CUfunction)tmp_kernel;
  {
    std::lock_guard<std::mutex> lock(cache_mutex);
    if (tuning.find(cf) != tuning.end()) {
      nvtxRangePop();
      return tuning.at(cf);
    }
  }

  debug("=== performing tuning for %s with default block "
        "%d,%d,%d, sub-block "
        "%d,%d,%d, and expected grid size (%d, %d, %d)",
        name, preferred.x, preferred.y, preferred.z, preferred_sub.x,
        preferred_sub.y, preferred_sub.z, expected_grid.x, expected_grid.y,
        expected_grid.z);

  int grid;

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

  int warps = threads / 32;

  debug("maximum warps = %d", warps);
  assert(warps >= 1);

  tuned_kernel result(preferred, preferred_sub);

  dim3 block(1, 1, 32);

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
    std::lock_guard<std::mutex> lock(cache_mutex);
    tuning[cf] = result;
    nvtxRangePop();
    return tuning[cf];
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
    nearest_square = max(1, (int)(floor(sqrtf(threads / block.z))));
    if ((small_dims[0] || small_dims[1]) &&
        _est_waste(expected_grid,
                   dim3(nearest_square, nearest_square, block.z)) > 0.2) {
      debug("remaining small dimensions would lead to excess wasted "
            "grid points, not trying a square");
    } else {
      next_block = dim3(
          block.x > 1 ? min(block.x, nearest_square) : nearest_square,
          block.y > 1 ? min(block.y, nearest_square) : nearest_square, block.z);
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
      tried.emplace(std::make_tuple(next_block.y, next_block.x, next_block.z));

      // find largest grid for which (kernel max occupancy threads) %
      // (x*y*z) ~= 0
      for (int offset = 0; offset < 3 && best_eff < 0.99999; offset++) {
        int xy = max_block / block.z;

        int x = (int)(floor(sqrtf(xy)));
        while (x > 1 && (xy % x > 0))
          x--;

        next_block = dim3(x, (xy / x), block.z);

        if (tried.find(std::make_tuple(next_block.x, next_block.y,
                                       next_block.z)) == tried.end()) {
          // next_block.x = threads / next_block.z / next_block.y;
          debug("next guess is %d,%d,%d", next_block.x, next_block.y,
                next_block.z);
          _check_kernel(next_block, preferred_sub, builder, next_valid, tmp,
                        max_block, next_regs, next_eff);

          debug(" (occupancy %.2f%%)", 100. * next_eff);
          if (compare_options(best_eff, std::get<0>(result), next_eff,
                              next_block)) {
            if (next_valid) {
              best_eff = next_eff;
              result = tuned_kernel(next_block, preferred_sub);
            }
          }
          tried.emplace(
              std::make_tuple(next_block.x, next_block.y, next_block.z));
          tried.emplace(
              std::make_tuple(next_block.y, next_block.x, next_block.z));
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
    std::lock_guard<std::mutex> lock(cache_mutex);
    tuning[cf] = result;
  }
  nvtxRangePop();
  return result;
}
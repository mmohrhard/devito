#pragma once

#include <cuda.h>
#include <devito/jitify.hpp>
#include <functional>
#include <map>
#include <stdio.h>
#include <tuple>
#include <utility>

#define TUNING_DEBUGGING

#ifdef TUNING_DEBUGGING
#define debug_printf(...) fprintf(stderr, __VA_ARGS__)
#else
#define debug_printf(...)                                                      \
  {}
#endif

#define CudaChecked(f) _cudaChecked((f), __FILE__, __LINE__)
#define STRINGIFY(x) _stringify(x)
#define _stringify(x) #x

#define MAX_CUDA_DEVICES 16
#define PER_DEVICE_TEMP_GET(T, NAME, SIZE)                                     \
  T *NAME = nullptr;                                                           \
  static T *NAME##_device[MAX_CUDA_DEVICES] = {0};                             \
  {                                                                            \
    int device = 0;                                                            \
    cudaGetDevice(&device);                                                    \
    assert(device >= 0 && device < MAX_CUDA_DEVICES);                          \
    if (NAME##_device[device] == nullptr) {                                    \
      CudaChecked(cudaMalloc((void **)&NAME##_device[device], (SIZE)));        \
      cudaMemset(NAME##_device[device], 1, (SIZE));                            \
    }                                                                          \
    NAME = NAME##_device[device];                                              \
  }

#define PER_DEVICE_TEMP_DESTROY(NAME)                                          \
  {                                                                            \
    int device = 0;                                                            \
    cudaGetDevice(&device);                                                    \
    assert(device >= 0 && device < MAX_CUDA_DEVICES);                          \
    CudaChecked(cudaFree(NAME##_device[device]));                              \
    NAME##_device[device] = nullptr;                                           \
  }

#define PER_DEVICE_ARRAY_TEMP_DECLARE(NAME, ARRAYTYPE)                         \
  static ARRAYTYPE *NAME##_device[MAX_CUDA_DEVICES] = {0};

#define PER_DEVICE_ARRAY_TEMP_GET(NAME, NBYTES)                                \
  _allocTempArray(&NAME##_device[_cudaGetCurrentDevice()], NBYTES)

#define PER_DEVICE_ARRAY_TEMP_DESTROY(NAME)                                    \
  {                                                                            \
    int device = _cudaGetCurrentDevice();                                      \
    _freeTempArray(NAME##_device[device]);                                     \
    NAME##_device[device] = nullptr;                                           \
  }

#define pow(x, y) powf(x, y)

inline void _cudaChecked(cudaError_t err, const char *file, int line,
                         const char *extra = nullptr) {
  if (err != cudaSuccess) {
    err = cudaGetLastError();
    fprintf(stderr, "!!! CUDA Error in operator: %s:%d %s\n", file, line,
            cudaGetErrorString(err));
    exit(1);
  }
}

#define CudaCheckLaunch(grid, block)                                           \
  _cudaCheckKernelLaunch(grid, block, __FILE__, __LINE__)
inline void _cudaCheckKernelLaunch(dim3 grid, dim3 block, const char *file,
                                   int line) {
  cudaError_t err = cudaPeekAtLastError();
  if (err != cudaSuccess) {
    fprintf(stderr, "!!! CUDA Error after kernel launch: %s:%d %s\n", file,
            line, cudaGetErrorString(err));
    exit(1);
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

template <typename T> inline T *_allocTempArray(T **array_ptr, size_t nbytes) {
  if (*array_ptr == nullptr) {
    CudaChecked(cudaMallocHost((void **)array_ptr, sizeof(T)));
    memset((void *)(*array_ptr), 0, sizeof(T));
  }

  T *array = *array_ptr;
  if (array->nbytes != nbytes) {
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
                        bool cond = true, cudaStream_t stream = nullptr) {
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
    if (stream != nullptr && stream != cudaStreamDefault)
      CudaChecked(cudaMemcpyAsync(dst, src, size, kind, stream));
    else {
      debug_printf("transferring %s..\n",
                   kind == cudaMemcpyHostToDevice ? "H->D" : "D->H");
      CudaChecked(cudaMemcpy(dst, src, size, kind));
    }
  }
}

template <typename T>
void prepareDataObject(T *obj, size_t size, bool copyIn = true,
                       cudaStream_t stream = nullptr) {
  if (!_cudaPtrIsManaged(obj->data)) {
    if (obj->device_data == nullptr) {
      CudaChecked(cudaMalloc((void **)&obj->device_data, size));
      obj->operator_allocated = 1;
    }
    transferDataObject(cudaMemcpyHostToDevice, obj, size, copyIn, stream);
  }
}

template <typename T> void destroyDataObject(T *obj, bool del = true) {
  if (del && obj->operator_allocated) {
    CudaChecked(cudaFree(obj->device_data));
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

static bool
_check_kernel(const dim3 &block, const dim3 &sub_block,
              std::function<jitify::KernelInstantiation(dim3, dim3)> &builder,
              bool &is_valid, float &est_occupancy) {
  // verify that it works by checking occupancy with the new block size
  auto k = builder(block, sub_block);

  int grid = 0;
  int max_block = 0;

  is_valid = false;
  est_occupancy = 0.f;

  CUresult res = cuOccupancyMaxPotentialBlockSize(&grid, &max_block,
                                                  (CUfunction)k, nullptr, 0, 0);
  if (res == 0) {
    debug_printf("rebuilt kernel is valid, has max occupancy at %d threads "
                 "(versus %d we calculated)\n",
                 max_block, block.x * block.y * block.z);
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

static tuned_kernel
performTuning(tuningDict &tuning, const char *name, dim3 preferred,
              dim3 preferred_sub, dim3 expected_grid, int max_block_dimension,
              std::function<jitify::KernelInstantiation(dim3, dim3)> builder) {

  // Calculate the maximum possible occupancy for the preferred block size
  auto tmp_kernel = builder(preferred, preferred_sub);
  CUfunction cf = (CUfunction)tmp_kernel;
  if (tuning.find(cf) != tuning.end()) {
    return tuning.at(cf);
  }

  debug_printf(
      "=== performing tuning for %s with default block %d,%d,%d, sub-block "
      "%d,%d,%d, and expected grid size (%d, %d, %d)\n",
      name, preferred.x, preferred.y, preferred.z, preferred_sub.x,
      preferred_sub.y, preferred_sub.z, expected_grid.x, expected_grid.y,
      expected_grid.z);

  int grid;
  int max_block;
  CUresult res =
      cuOccupancyMaxPotentialBlockSize(&grid, &max_block, cf, nullptr, 0, 0);

  if (res != 0) {
    fprintf(stderr,
            "!!! invalid kernel detected! could not calculate occupancy for %s "
            "with blocksize (%d, %d, %d), sub-block size (%d, %d, %d)!\n",
            name, preferred.x, preferred.y, preferred.z, preferred_sub.x,
            preferred_sub.y, preferred_sub.z);
    assert(false);
  }

  debug_printf("maximum occupancy is at %d threads with minimum %d grid size\n",
               max_block, grid);

  int threads = max_block;

  int last_dim = dim3_get(expected_grid, max_block_dimension - 1);

  int warps = threads / 32;

  debug_printf("maximum warps = %d\n", warps);
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
    debug_printf("expected grid (%d, %d, %d) has one or more small dimensions "
                 "(<32 or >75%% expected waste)\n",
                 expected_grid.x, expected_grid.y, expected_grid.z);

  // dim3 sub(1, 1, 1);
  int nearest_square = 0;
  // int diff = threads - (preferred.x * preferred.y * preferred.z);
  // int best_diff = diff;

  // Assume the initial input is 'best' until we know otherwise
  float best_eff = _est_efficiency(max_block, preferred);
  debug_printf("base efficiency is %.2f\n", best_eff);
  float next_eff = 0.f;
  bool next_valid = false;
  float tmp = 0.f;

  dim3 next_block;

  switch (max_block_dimension) {
  case 1:
    // We're done, 1D is a simple case
    debug_printf("kernel is 1D and maximum block size is %d\n", max_block);
    std::get<0>(result).x = max_block;
    break;

  case 2:
    // for now, just use preferred size, even though that's suboptimal nearly
    // always

    break;

  case 3:
    // Force small dimensions to have small block sizes, even if it means warp
    // divergence
    for (int d = 0; d < 3; d++) {
      if (small_dims[d])
        dim3_set(block, d, 4);
    }

    // try and make a better guess
    nearest_square = max(1, (int)(floor(sqrtf(threads / block.z))));
    if ((small_dims[0] || small_dims[1]) &&
        _est_waste(expected_grid,
                   dim3(nearest_square, nearest_square, block.z)) > 0.2) {
      debug_printf("remaining small dimensions would lead to excess wasted "
                   "grid points, not trying a square\n");
    } else {
      next_block = dim3(
          block.x > 1 ? min(block.x, nearest_square) : nearest_square,
          block.y > 1 ? min(block.y, nearest_square) : nearest_square, 
          block.z);
      next_eff = _est_efficiency(threads, next_block);

      _check_kernel(next_block, preferred_sub, builder, next_valid, tmp);
      if (next_valid) {
        debug_printf(
            "nearest square attempt would be %d,%d,%d (efficiency %.2f%%)\n",
            next_block.x, next_block.y, next_block.z, 100. * next_eff);

        if (next_eff > best_eff) {
          result = tuned_kernel(next_block, preferred_sub);
          best_eff = next_eff;
        }
      }

      next_block = dim3(1, block.y > 1 ? min(block.y, nearest_square & ~1) : (nearest_square & ~1), block.z);
      next_block.x = threads / next_block.z / next_block.y;
      next_eff = _est_efficiency(threads, next_block);
      if (next_eff > best_eff) {
        debug_printf("next guess is %d,%d,%d (efficiency %.2f%%)\n",
                     next_block.x, next_block.y, next_block.z, 100. * next_eff);
        _check_kernel(next_block, preferred_sub, builder, next_valid, tmp);
        if (next_valid) {
          best_eff = next_eff;
          result = tuned_kernel(next_block, preferred_sub);
        }
      }
    }

    break;
  default:
    assert(false);
  }

  debug_printf("selected block (%d, %d, %d), sub_block (%d, %d, %d)\n",
               std::get<0>(result).x, std::get<0>(result).y,
               std::get<0>(result).z, std::get<1>(result).x,
               std::get<1>(result).y, std::get<1>(result).z);
  tuning[cf] = result;
  return tuning[cf];
}
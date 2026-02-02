#ifndef _DEVITO_CUDA_KERNELS_H
#define _DEVITO_CUDA_KERNELS_H

#include <cuda_runtime.h>
#include <cstddef>
#include <map>
#include <memory>
#include <set>

#include <devito/types.hpp>
#include <devito/jitify.hpp>
#include <devito/logging.hpp>
#include <devito/errors.hpp>
#include <devito/util.hpp>


using tuned_kernel = std::pair<dim3, dim3>;

typedef std::map<std::shared_ptr<jitify::detail::CUDAKernel>, tuned_kernel>
    tuningDict;


namespace devito {
namespace cuda {

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

  // Convert our minimum values to block-aligned minimums
  dim3 block_min = dim3((mins.x / block_sizes.x) * block_sizes.x,
                        (mins.y / block_sizes.y) * block_sizes.y,
                        (mins.z / block_sizes.z) * block_sizes.z);

  // Resize the grid based on the block sizes and offsets
  dim3 grid(ceil_div(maxs.x - block_min.x, block_sizes.x),
            ceil_div(maxs.y - block_min.y, block_sizes.y),
            ceil_div(maxs.z - block_min.z, block_sizes.z));

  if (grid.x >= 1 && grid.y >= 1 && grid.z >= 1) {
    if (!_cudaChecked(
            (cudaError_t)(kernel.configure(grid, tb, 0, stream)
                              .launch(block_min, std::forward<Args>(args)...)),
            file, line)) {
      return -1;
    }
  }

  return 0;
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

} // namespace cuda
} // namespace devito

#endif // _DEVITO_CUDA_KERNELS_H

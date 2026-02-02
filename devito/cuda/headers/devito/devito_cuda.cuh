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

#include <devito/errors.hpp>
#include <devito/kernels.hpp>
#include <devito/logging.hpp>
#include <devito/memory.hpp>
#include <devito/mpi.hpp>
#include <devito/profiling.hpp>
#include <devito/types.hpp>

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
          cudaMemsetAsync(NAME##_device[device], 0, (SIZE), kernel_stream));   \
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
    if (devito::cuda::_allocTempArray<DTYPE>(                                  \
            &NAME##_device[devito::cuda::_cudaGetCurrentDevice()],             \
            STRINGIFY(NAME), {__VA_ARGS__}) != 0)                              \
      return -1;                                                               \
    NAME = NAME##_device[devito::cuda::_cudaGetCurrentDevice()];               \
  }

#define PER_DEVICE_ARRAY_TEMP_DESTROY(NAME)                                    \
  {                                                                            \
    int device = devito::cuda::_cudaGetCurrentDevice();                        \
    devito::cuda::_freeTempArray(NAME##_device[device]);                       \
    NAME##_device[device] = nullptr;                                           \
  }

#define pow(x, y) powf(x, y)

/**
 * Miscellaneous helpers
 */

#define ENSURE_STREAM_PRIO(NAME, PRIORITY)                                     \
  static cudaStream_t NAME##_devs[MAX_CUDA_DEVICES] = {0};                     \
  (void)NAME;                                                                  \
  {                                                                            \
    int device = devito::cuda::_cudaGetCurrentDevice();                        \
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
    int cur = devito::cuda::_cudaGetCurrentDevice();                           \
    if (cur != dev && dev != -1)                                               \
      warn("!!! CUDA issue: expected device %d, current device %d", dev, cur); \
    cudaSetDevice(dev);                                                        \
  }

#define prepareDataObject(NAME, ...)                                           \
  devito::cuda::_prepareDataObject(                                            \
      NAME, STRINGIFY(NAME), "prepareDataObject(" #NAME ")", __VA_ARGS__)

#define destroyDataObject(NAME, ...)                                           \
  {                                                                            \
    if (devito::cuda::_destroyDataObject(NAME, STRINGIFY(NAME), __VA_ARGS__) < \
        0)                                                                     \
      return -1;                                                               \
  }

#define transferDataObject(...)                                                \
  {                                                                            \
    if (devito::cuda::transferDataObject(__VA_ARGS__) < 0)                     \
      return -1;                                                               \
  }

#define launchKernel(KNAME, ...)                                               \
  {                                                                            \
    if (devito::cuda::_launchKernel(__FILE__, __LINE__, STRINGIFY(KNAME),      \
                                    KNAME##_tuned, KNAME##_tune,               \
                                    __VA_ARGS__) < 0) {                        \
      return -1;                                                               \
    }                                                                          \
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

#define DEVITO_CUDA_PROLOGUE()                                                 \
  devito::cuda::CudaSectionTimers _cuda_section_timers(kernel_stream);         \
  nvtxRangePush(__FUNCTION__);

#define DEVITO_CUDA_EPILOGUE()                                                 \
  do {                                                                         \
    _cuda_section_timers.resolveTimers();                                      \
    if (devicerm || updatehost) {                                              \
      CudaChecked(cudaStreamSynchronize(kernel_stream));                       \
      nvtxRangePop();                                                          \
    }                                                                          \
  } while (0);

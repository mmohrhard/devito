#ifndef _DEVITO_CUDA_ERRORS_H
#define _DEVITO_CUDA_ERRORS_H

#include <cuda_runtime.h>
#include <nccl.h>
#include <string>

#ifndef OPERATOR_STANDALONE
#include <Python.h>
#endif

#include <devito/logging.hpp>

namespace devito {
namespace cuda {

/**
 * @brief Raise a Python error.
 *
 * @param format The format string
 * @param args Format string arguments
 */
template <typename... Args>
__host__ void acquire_gil_and_raise_error(const std::string &format, Args... args) {
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
__host__ inline bool _cudaChecked(cudaError_t err, const char *file, int line,
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

__host__ inline bool _ncclChecked(ncclResult_t err, const char *file, int line,
                         const char *extra = nullptr) {
  // We may need this if we switch to fully-async NCCL, but I'm not sure
  // what that actually buys us?
  //
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
__host__ inline void _cudaCheckKernelLaunch(dim3 grid, dim3 block, const char *file,
                                   int line) {
  cudaError_t err = cudaPeekAtLastError();
  if (err != cudaSuccess) {
    critical("!!! CUDA Error after kernel launch: %s:%d %s", file, line,
             cudaGetErrorString(err));
  }
}

} // namespace cuda
} // namespace devito

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
    if (!devito::cuda::_cudaChecked((cudaError_t)(f), __FILE__, __LINE__))     \
      return -1;                                                               \
  }
#define CudaCheckedEx(f, msg)                                                  \
  {                                                                            \
    if (!devito::cuda::_cudaChecked((cudaError_t)(f), __FILE__, __LINE__,      \
                                    msg))                                      \
      return -1;                                                               \
  }
#define NcclChecked(f)                                                         \
  {                                                                            \
    if (!devito::cuda::_ncclChecked((ncclResult_t)(f), __FILE__, __LINE__))    \
      return -1;                                                               \
  }
#define NcclCheckedEx(f, msg)                                                  \
  {                                                                            \
    if (!devito::cuda::_ncclChecked((ncclResult_t)(f), __FILE__, __LINE__,     \
                                    msg))                                      \
      return -1;                                                               \
  }

#endif // _DEVITO_CUDA_ERRORS_H

#ifndef _DEVITO_CUDA_MPI_H
#define _DEVITO_CUDA_MPI_H

// MPI bits bloew
#ifdef MPI_VERSION

#include <cstddef>
#include <utility>
#include <cuda_runtime.h>
#include <nccl.h>

#include <devito/errors.hpp>
#include <devito/logging.hpp>

namespace devito {
namespace cuda {

// Multi-function halo update using NCCL
//
// Basically just exists to avoid Devito spitting out a ton of nearly
// identical halo exchange functions.
//
// That, and we can do C++-y things more easily in a C++ header than in
// Devito's AST/IR format.
//
template <typename TDataobj, typename TMPIMsg>
static int async_multi_haloupdate(
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
#ifdef DEVITO_CUDA_VERBOSE_GATHER_SCATTER
        debug("devito_cuda_async_multi_haloupdate: func=%p, buf=%p, "
              "sizes=(%d,%d,%d), "
              "otime=%d, ofs=(%d,%d,%d), stream=%p, nsizes=%d",
              function, msg[i].bufg, msg[i].sizes[0], msg[i].sizes[1],
              msg[i].sizes[2], otime, msg[i].ofsg[0], msg[i].ofsg[1],
              msg[i].ofsg[2], (void *)kernel_stream, msg[i].nsizes);
#endif
        Checked(async_gather_4d<TDataobj>(
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
async_multi_halowait(std::initializer_list<TDataobj> functions,
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
#ifdef DEVITO_CUDA_VERBOSE_GATHER_SCATTER
        debug("devito_cuda_async_multi_halowait: func=%p, buf=%p, "
              "sizes=(%d,%d,%d), "
              "otime=%d, ofs=(%d,%d,%d), stream=%p, nsizes=%d",
              function, msg[i].bufs, msg[i].sizes[0], msg[i].sizes[1],
              msg[i].sizes[2], otime, msg[i].ofsg[0], msg[i].ofsg[1],
              msg[i].ofsg[2], (void *)kernel_stream, msg[i].nsizes);
#endif
        Checked(async_scatter_4d<TDataobj>(
            function, (float *)msg[i].bufs, msg[i].sizes[0], msg[i].sizes[1],
            msg[i].sizes[2], otime, msg[i].ofsg[0], msg[i].ofsg[1],
            msg[i].ofsg[2], kernel_stream));
      }
    }
  }

  return 0;
}

} // namespace cuda
} // namespace devito

#endif // MPI_VERSION
#endif // _DEVITO_CUDA_MPI_H

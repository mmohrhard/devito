#ifndef _DEVITO_CUDA_MPI_H
#define _DEVITO_CUDA_MPI_H

// MPI bits bloew

#include <array>
#include <cstddef>
#include <cuda_runtime.h>
#include <nccl.h>
#include <utility>

#include <devito/errors.hpp>
#include <devito/logging.hpp>
#include <devito/memory.hpp>
#include <devito/types.hpp>

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
template <typename TDataobj, typename TMPIMsg, int ncomms>
static int async_multi_haloupdate(std::initializer_list<TDataobj> functions,
                                  std::initializer_list<TMPIMsg> msgs,
                                  int otime, cudaStream_t nccl_stream,
                                  cudaStream_t kernel_stream,
                                  ncclComm_t nccl_comm) {

  assert(functions.size() == msgs.size());
  cudaEvent_t update_ev = nullptr;
  CudaChecked(cudaEventCreateWithFlags(&update_ev, cudaEventDisableTiming));

  TDataobj *fptr = const_cast<TDataobj *>(functions.begin());
  TMPIMsg *mptr = const_cast<TMPIMsg *>(msgs.begin());

  // issue gathers on the kernel stream
  for (int f = 0; f < functions.size(); f++) {
    TDataobj function = fptr[f];
    TMPIMsg msg = mptr[f];
    for (int i = 0; i < ncomms; i++) {
      if (msg[i].torank != MPI_PROC_NULL) {
#ifdef DEVITO_CUDA_VERBOSE_GATHER_SCATTER
        debug("devito_cuda_async_multi_haloupdate: func=%s(%p), ncomm=%d, "
              "buf=%p, "
              "size=(%d,%d,%d) "
              "extent=(%d,%d,%d), "
              "otime=%d, ofs=(%d,%d,%d), stream=%p, nsizes=%d",
              function->name, function, i, msg[i].bufg, function->size[1],
              function->size[2], function->size[3], msg[i].sizes[0],
              msg[i].sizes[1], msg[i].sizes[2], otime, msg[i].ofsg[0],
              msg[i].ofsg[1], msg[i].ofsg[2], (void *)kernel_stream,
              msg[i].nsizes);
#endif
        Checked(async_gather_4d<TDataobj>(
            (float *)msg[i].bufg, function, msg[i].sizes[0], msg[i].sizes[1],
            msg[i].sizes[2], otime, msg[i].ofsg[0], msg[i].ofsg[1],
            msg[i].ofsg[2], kernel_stream));
      }
    }
  }

  // pivot to the nccl stream for the actual device exchanges
  CudaChecked(cudaEventRecord(update_ev, kernel_stream));
  CudaChecked(cudaStreamWaitEvent(nccl_stream, update_ev));
  CudaChecked(cudaEventDestroy(update_ev));

  NcclChecked(ncclGroupStart());

  for (int f = 0; f < functions.size(); f++) {
    TDataobj function = fptr[f];
    TMPIMsg msg = mptr[f];
    for (int i = 0; i < ncomms; i++) {
      // NCCL isn't MPI, so we omit send/recv to MPI_PROC_NULL ranks explicitly
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

template <typename TDataobj, typename TMPIMsg, int ncomms>
static int async_multi_halowait(std::initializer_list<TDataobj> functions,
                                std::initializer_list<TMPIMsg> msgs, int otime,
                                cudaStream_t nccl_stream,
                                cudaStream_t kernel_stream,
                                ncclComm_t nccl_comm) {
  assert(functions.size() == msgs.size());
  cudaEvent_t update_ev = nullptr;
  CudaChecked(cudaEventCreateWithFlags(&update_ev, cudaEventDisableTiming));

  // Pivot back to the kernel stream from the halo exchanges on the nccl stream
  CudaChecked(cudaEventRecord(update_ev, nccl_stream));
  CudaChecked(cudaStreamWaitEvent(kernel_stream, update_ev));
  CudaChecked(cudaEventDestroy(update_ev));

  TDataobj *fptr = const_cast<TDataobj *>(functions.begin());
  TMPIMsg *mptr = const_cast<TMPIMsg *>(msgs.begin());

  // issue the halo scatters on the kernel stream
  for (int f = 0; f < functions.size(); f++) {
    TDataobj function = fptr[f];
    TMPIMsg msg = mptr[f];
    for (int i = 0; i < ncomms; i++) {
      if (msg[i].fromrank != MPI_PROC_NULL) {
#ifdef DEVITO_CUDA_VERBOSE_GATHER_SCATTER
        debug("devito_cuda_async_multi_halowait: func=%s(%p), buf=%p, "
              "size=(%d,%d,%d) "
              "extent=(%d,%d,%d), "
              "otime=%d, ofs=(%d,%d,%d), stream=%p, nsizes=%d",
              function->name, function, msg[i].bufs, function->size[1],
              function->size[2], function->size[3], msg[i].sizes[0],
              msg[i].sizes[1], msg[i].sizes[2], otime, msg[i].ofss[0],
              msg[i].ofss[1], msg[i].ofss[2], (void *)nccl_stream,
              msg[i].nsizes);
#endif
        Checked(async_scatter_4d<TDataobj>(
            function, (float *)msg[i].bufs, msg[i].sizes[0], msg[i].sizes[1],
            msg[i].sizes[2], otime, msg[i].ofss[0], msg[i].ofss[1],
            msg[i].ofss[2], kernel_stream));
      }
    }
  }

  return 0;
}

} // namespace cuda
} // namespace devito

#endif // _DEVITO_CUDA_MPI_H

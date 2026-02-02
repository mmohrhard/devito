#ifndef _DEVITO_CUDA_MEMORY_H
#define _DEVITO_CUDA_MEMORY_H

#include "nvtx3/nvToolsExt.h"
#include <cstddef>
#include <cuda_runtime.h>

#include <devito/types.hpp>
#include <devito/errors.hpp>
#include <devito/logging.hpp>
#include <devito/util.hpp>

namespace devito {
namespace cuda {

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

  if (size == 0) {
    return 0; // Nothing to do
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
      if (size > 0) {
        debug("allocating %llu bytes for %s on device %d", size, name, device);
        CudaChecked(cudaMallocAsync((void **)&obj->device_data, size, stream));
        obj->operator_allocated = 1;
      }
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
    if (size > 0) {
      ret = transferDataObject(cudaMemcpyHostToDevice, obj, size, copyIn,
                               stream, name);
    }
  }
  nvtxRangePop();
  return ret;
}

template <typename T>
int _destroyDataObject(T *obj, const char *name, bool del = true,
                       cudaStream_t stream = nullptr) {
  if (del && obj->operator_allocated && obj->device_data != nullptr) {
    debug("freeing %s", name);
    CudaChecked(cudaFreeAsync(obj->device_data, stream));
    obj->device_data = nullptr;
  }

  return 0;
}

/// Get a human-readable string describing the owner of a CUDA pointer
template <typename T> inline static char *_cudaPtrOwnerString(T *ptr) {
  struct cudaPointerAttributes attr = cudaPointerAttributes{};
  int ret = cudaPointerGetAttributes(&attr, (const void *)ptr);
  if (ret != cudaSuccess) {
    return (char *)"Failed";
  } else if (attr.type == cudaMemoryType::cudaMemoryTypeDevice) {
    return (char *)"Device";
  } else if (attr.type == cudaMemoryType::cudaMemoryTypeHost) {
    return (char *)"Host";
  } else if (attr.type == cudaMemoryType::cudaMemoryTypeManaged) {
    return (char *)"Managed";
  } else if (attr.type == cudaMemoryType::cudaMemoryTypeUnregistered) {
    return (char *)"Unregistered Host Memory";
  } else {
    return (char *)"Unknown";
  }
}

/// Perform an asynchronous gather from a 4D data object into a contiguous
/// buffer for MPI/NCCL interchange purposes
template <typename TDataobj>
static int async_gather_4d(void *__restrict buf, TDataobj dataobj,
                           const int x_sz, const int y_sz, const int z_sz,
                           const int w_ofs, const int x_ofs, const int y_ofs,
                           const int z_ofs, cudaStream_t stream) {
  struct cudaMemcpy3DParms copy_params = {0};
  if (!_cudaPtrIsDevicePtr(dataobj->device_data)) {
    critical(
        "!!! devito::cuda::async_gather_4d: dataobj->device_data %p is not "
        "a device pointer (is %s)",
        dataobj->device_data, _cudaPtrOwnerString(dataobj->device_data));
    return -1;
  }
  if (!_cudaPtrIsDevicePtr(buf)) {
    debug(
        "!!! devito::cuda::async_gather_4d: buf %p is not a device pointer (is "
        "%s)",
        buf, _cudaPtrOwnerString(buf));
    return -1;
  }

  assert(dataobj->rank == 4);

  int rank = dataobj->rank;

  size_t pitch = dataobj->size[rank - 1] * dataobj->element_size;
  size_t row_size = dataobj->size[rank - 1] * dataobj->element_size;
  size_t ptr_offset = dataobj->element_size * w_ofs * dataobj->size[rank - 3] *
                      dataobj->size[rank - 2] * dataobj->size[rank - 1];

#ifdef DEVITO_CUDA_VERBOSE_GATHER_SCATTER
  debug("devito::cuda::async_gather_4d: buf=%p, dataobj=%p, x_sz=%d, y_sz=%d, "
        "z_sz=%d (%d bytes), "
        "w_ofs=%d, "
        "x_ofs=%d, y_ofs=%d, z_ofs=%d, stream=%p, element_size=%d",
        buf, dataobj->device_data, x_sz, y_sz, z_sz,
        dataobj->element_size * z_sz, w_ofs, x_ofs, y_ofs, z_ofs,
        (void *)stream, dataobj->element_size);
  debug("devito::cuda::async_gather_4d: src pitch=%llu, row_size=%llu, "
        "ptr_offset=0x%llx",
        pitch, row_size, ptr_offset);
  debug("devito::cuda::async_gather_4d: src ptr alignment=%llu",
        1 << __builtin_ctzll(
            (uintptr_t)(((char *)(dataobj->device_data)) + ptr_offset)));
  debug("devito::cuda::async_gather_4d: dst ptr alignment=%llu",
        1 << __builtin_ctzll((uintptr_t)(buf)));
  debug("devito::cuda::async_gather_4d: dst pitch=%llu, row_size=%llu",
        dataobj->element_size * (z_sz), dataobj->element_size * (z_sz));
#endif

  copy_params.srcPtr = make_cudaPitchedPtr(
      reinterpret_cast<char *>(dataobj->device_data) + ptr_offset, pitch,
      row_size, dataobj->size[rank - 2]);
  copy_params.srcPos =
      make_cudaPos(dataobj->element_size * z_ofs, y_ofs, x_ofs);
  copy_params.dstPtr = make_cudaPitchedPtr(buf, dataobj->element_size * z_sz,
                                           dataobj->element_size * z_sz, y_sz);
  copy_params.dstPos = make_cudaPos(0, 0, 0);
  copy_params.extent =
      make_cudaExtent(dataobj->element_size * z_sz, y_sz, x_sz);
  copy_params.kind = cudaMemcpyDeviceToDevice;
  CudaChecked(cudaMemcpy3DAsync(&copy_params, stream));
  return cudaSuccess;
}

/// Perform an asynchronous scatter from a contiguous buffer into a 4D data
/// object
template <typename TDataobj>
static int async_scatter_4d(TDataobj dataobj, const void *const __restrict buf,
                            const int x_sz, const int y_sz, const int z_sz,
                            const int w_ofs, const int x_ofs, const int y_ofs,
                            const int z_ofs, cudaStream_t stream) {
  struct cudaMemcpy3DParms copy_params = {0};
  if (!_cudaPtrIsDevicePtr(dataobj->device_data)) {
    critical(
        "!!! devito::cuda::async_scatter_4d: dataobj->device_data %p is not "
        "a device pointer (is %s)",
        dataobj->device_data, _cudaPtrOwnerString(dataobj->device_data));
    return -1;
  }
  if (!_cudaPtrIsDevicePtr(buf)) {
    debug("!!! devito::cuda::async_scatter_4d: buf %p is not a device pointer "
          "(is %s)",
          buf, _cudaPtrOwnerString(buf));
    return -1;
  }

  assert(dataobj->rank == 4);
  int rank = dataobj->rank;
  size_t pitch = dataobj->size[rank - 1] * dataobj->element_size;
  size_t row_size = dataobj->size[rank - 1] * dataobj->element_size;
  size_t ptr_offset = dataobj->element_size * w_ofs * dataobj->size[rank - 3] *
                      dataobj->size[rank - 2] * dataobj->size[rank - 1];

#ifdef DEVITO_CUDA_VERBOSE_GATHER_SCATTER
  debug("devito::cuda::async_scatter_4d: buf=%p, dataobj=%p, x_sz=%d, y_sz=%d, "
        "z_sz=%d (%d bytes), "
        "w_ofs=%d, "
        "x_ofs=%d, y_ofs=%d, z_ofs=%d, stream=%p, element_size=%d",
        buf, dataobj->device_data, x_sz, y_sz, z_sz,
        dataobj->element_size * z_sz, w_ofs, x_ofs, y_ofs, z_ofs,
        (void *)stream, dataobj->element_size);
  debug("devito::cuda::async_scatter_4d: src pitch=%llu, row_size=%llu, "
        "ptr_offset=0x%llx",
        pitch, row_size, ptr_offset);
  debug("devito::cuda::async_scatter_4d: src ptr alignment=%llu",
        1 << __builtin_ctzll(
            (uintptr_t)(((char *)(dataobj->device_data)) + ptr_offset)));
  debug("devito::cuda::async_scatter_4d: dst ptr alignment=%llu",
        1 << __builtin_ctzll((uintptr_t)(buf)));
  debug("devito::cuda::async_scatter_4d: dst pitch=%llu, row_size=%llu",
        dataobj->element_size * z_sz, dataobj->element_size * z_sz);
#endif
  copy_params.dstPtr = make_cudaPitchedPtr(
      reinterpret_cast<char *>(dataobj->device_data) + ptr_offset, pitch,
      row_size, dataobj->size[rank - 2]);
  copy_params.dstPos =
      make_cudaPos(dataobj->element_size * z_ofs, y_ofs, x_ofs);
  copy_params.srcPtr =
      make_cudaPitchedPtr(const_cast<void *>(buf), dataobj->element_size * z_sz,
                          dataobj->element_size * z_sz, y_sz);
  copy_params.srcPos = make_cudaPos(0, 0, 0);
  copy_params.extent =
      make_cudaExtent(dataobj->element_size * z_sz, y_sz, x_sz);
  copy_params.kind = cudaMemcpyDeviceToDevice;

  CudaChecked(cudaMemcpy3DAsync(&copy_params, stream));
  return cudaSuccess;
}



/// Asynchronous savebuffer copy from device to host or vice versa
///
/// Eventually, we'll de inline compression/decompression in here
/// as required.
template <typename TSrc, typename TDst>
static int async_buffer_copy_slice(TDst *dst, const TSrc *src, const int t_dst,
                                   const int t_src, cudaMemcpyKind kind,
                                   cudaStream_t stream, dim4 src_offset = dim4(), dim4 dst_offset = dim4(), dim4 extent = dim4(0,0,0,0)) {
  // Copy a single time slice of 3D data from src to dst
  struct cudaMemcpy3DParms copy_params = {0};
  assert(src->rank == 4 || t_src == 0);
  assert(dst->rank == 4 || t_dst == 0);
  int srank = src->rank;
  int drank = dst->rank;


  size_t pitch_src = src->size[srank - 1] * src->element_size;
  size_t pitch_dst = dst->size[drank - 1] * dst->element_size;

  // For Devito, the row size is also the pitch - we just
  // copy all the padding to make it easier
  size_t row_size_src = src->size[srank - 1] * src->element_size;
  size_t row_size_dst = dst->size[drank - 1] * dst->element_size;

  // Calculate pointer offsets in the fourth dimension
  // since CUDA doesn't directly handle anything beyond the
  // third dimension
  size_t ptr_offset_src = src->element_size * t_src * src->size[srank - 3] *
                          src->size[srank - 2] * src->size[srank - 1];
  size_t ptr_offset_dst = dst->element_size * t_dst * dst->size[drank - 3] *
                          dst->size[drank - 2] * dst->size[drank - 1];

  char *src_ptr = reinterpret_cast<char *>(select_ptr(src, kind == cudaMemcpyDeviceToHost));
  char *dst_ptr = reinterpret_cast<char *>(select_ptr(dst, kind == cudaMemcpyHostToDevice));

  copy_params.srcPtr = make_cudaPitchedPtr(src_ptr + ptr_offset_src, pitch_src,
                                           row_size_src, src->size[srank - 2]);
  copy_params.srcPos = make_cudaPos(0, 0, 0);

  copy_params.dstPtr = make_cudaPitchedPtr(dst_ptr + ptr_offset_dst, pitch_dst,
                                           row_size_dst, dst->size[drank - 2]);
  copy_params.dstPos = make_cudaPos(0, 0, 0);

  // Copy the entire 3D slice
  copy_params.extent =
      make_cudaExtent(src->element_size * (src->size[srank - 1]),
                      src->size[srank - 2], src->size[srank - 3]);
  // This may be used either for inbound or outbound copies
  copy_params.kind = kind;
  CudaChecked(cudaMemcpy3DAsync(&copy_params, stream));
  return cudaSuccess;
}

namespace pair_iterators {
template <typename T1, typename T2>
std::pair<T1, T2> operator++(std::pair<T1, T2> &it) {
  ++it.first;
  ++it.second;
  return it;
}
} // namespace pair_iterators

/* ========================================================================= */
/*
 * Bulk data object initialization / transfer / teardown
 */

template <typename TDataobj>
static int async_d2h_destroy_many(std::initializer_list<TDataobj> dataobjs,
                                  bool devicerm, bool updatehost,
                                  cudaStream_t stream) {
  nvtxRangePush("devito::cuda::async_d2h_destroy_many");
  int ret = cudaSuccess;
  for (auto &dobj : dataobjs) {
    char *name = dobj->name != nullptr ? dobj->name : "unnamed";
    char nvtxRange[256];
    memset(nvtxRange, 0, 256);
    snprintf(nvtxRange, 256, "transferDataObj(%s)", name);
    nvtxRangePush(nvtxRange);
    if (transferDataObject(cudaMemcpyDeviceToHost, dobj, 0, updatehost, stream,
                           name) < 0 ||
        _destroyDataObject(dobj, name, devicerm, stream) < 0) {
      ret = -1;
      nvtxRangePop();
      break;
    }
    nvtxRangePop();
  }
  nvtxRangePop();
  return ret;
}

template <typename TDataobj>
static int async_h2d_prepare_many(std::initializer_list<TDataobj> dataobjs,
                                  bool devicecreate, bool updatedevice,
                                  cudaStream_t stream) {
  nvtxRangePush("devito::cuda::async_h2d_prepare_many");
  char nvtxRange[256];
  int ret = 0;
  for (auto &dobj : dataobjs) {
    char *name = dobj->name != nullptr ? dobj->name : "unnamed";
    memset(nvtxRange, 0, 256);
    snprintf(nvtxRange, 256, "prepareDataObject(%s)", name);
    nvtxRangePush(nvtxRange);
    if (_prepareDataObject(dobj, name, nvtxRange, 0,
                           devicecreate || updatedevice, stream) < 0) {
      ret = -1;
      nvtxRangePop();
      break;
    }
    nvtxRangePop();
  }
  nvtxRangePop();
  return 0;
}

} // namespace cuda
} // namespace devito

#endif // _DEVITO_CUDA_MEMORY_H

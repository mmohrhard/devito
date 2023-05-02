#define MAX_CUDA_DEVICES 16
#define PER_DEVICE_TEMP_GET(T, NAME, SIZE)                                     \
  static T *NAME##_device[MAX_CUDA_DEVICES] = {0};                             \
  {                                                                            \
    int device = 0;                                                            \
    cudaGetDevice(&device);                                                    \
    assert(device >= 0 && device < MAX_CUDA_DEVICES);                          \
    if (NAME##_device[device] == nullptr) {                                    \
      CudaChecked(cudaMalloc((void **)&NAME##_device[device], (SIZE)));        \
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

#define PER_DEVICE_ARRAY_TEMP_GET(NAME, NBYTES)                     \
  _allocTempArray(&NAME##_device[_cudaGetCurrentDevice()], NBYTES)

#define PER_DEVICE_ARRAY_TEMP_DESTROY(NAME)                                    \
  {                                                                            \
    int device = _cudaGetCurrentDevice();                                      \
    _freeTempArray(NAME##_device[device]);                                     \
    NAME##_device[device] = nullptr;                                           \
  }

#define CudaChecked(f) _cudaChecked((f), __FILE__, __LINE__)
#define STRINGIFY(x) _stringify(x)
#define _stringify(x) #x

void _cudaChecked(cudaError_t err, const char *file, int line,
                  const char *extra = nullptr) {
  if (err != cudaSuccess) {
    err = cudaGetLastError();
    fprintf(stderr, "!!! CUDA Error in operator: %s:%d %s\n", file, line,
            cudaGetErrorString(err));
    exit(1);
  }
}

#define CudaCheckLaunch(f) _cudaCheckKernelLaunch((f), )
void _cudaCheckKernelLaunch(dim3 grid, dim3 block, const char *file, int line) {
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

template <typename T> T *_allocTempArray(T **array_ptr, size_t nbytes) {
  if (*array_ptr == nullptr) {
    CudaChecked(cudaMallocHost((void **)array_ptr, sizeof(T)));
    memset((void *)(*array_ptr), 0, sizeof(T));
  }

  T *array = *array_ptr;
  if (array->nbytes != nbytes) {
    _freeTempArrayData(array);
    CudaChecked(cudaMallocHost((void **)(&array->data), nbytes));
    CudaChecked(cudaMalloc((void **)(&array->device_data), nbytes));
    array->nbytes = nbytes;
  }

  return array;
}

uint64_t next_pow2(uint64_t x) {
  return (__builtin_popcount(x) == 1 || x == 1)
             ? x
             : 1 << (64 - __builtin_clzl(x - 1));
}

#define setupGrid(GRID, THREAD, X, Y, Z)                                       \
  _setupGrid(#GRID, GRID, THREAD, X, Y, Z)

// Heuristic thread block sizing based on the grid
// Not perfect, but it'll do for now
void _setupGrid(const char *gridName, dim3 &grid, dim3 &threadBlock,
                       int x_size, int y_size, int z_size) {
  //  long gp = x_size * y_size * z_size;
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

  grid = dim3((int)ceil((float)x_size / (float)threadBlock.x),
              (int)ceil((float)y_size / (float)threadBlock.y),
              (int)ceil((float)z_size / (float)threadBlock.z));
}

inline int _cudaGetCurrentDevice() {
  int device = -1;
  CudaChecked(cudaGetDevice(&device));
  return device;
}

template<typename T>
bool _cudaPointerIsAccessible(T* ptr) {
  struct cudaPointerAttributes attr = cudaPointerAttributes{};
  CudaChecked(cudaPointerGetAttributes(&attr, (const void *)ptr));
  return attr.devicePointer != NULL;
}

template <typename T> bool _cudaPtrIsManaged(T *ptr) {
    struct cudaPointerAttributes attr = cudaPointerAttributes{};
    CudaChecked(cudaPointerGetAttributes(&attr, (const void *)ptr));
    return attr.type == cudaMemoryType::cudaMemoryTypeManaged;
}

template <typename T> bool _cudaPtrIsDeviceAccessible(T* ptr) {
    if (ptr == nullptr)
        return false;

    struct cudaPointerAttributes attr = cudaPointerAttributes{};
    auto ret = cudaPointerGetAttributes(&attr, (const void *)ptr);
    return ret == cudaSuccess && attr.devicePointer != NULL;
}
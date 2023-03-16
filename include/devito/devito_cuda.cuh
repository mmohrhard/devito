#define CudaChecked(f) _cudaChecked((f), __FILE__, __LINE__)
#define STRINGIFY(x) _stringify(x)
#define _stringify(x) #x

void _cudaChecked(cudaError_t err, const char *file, int line,
                  const char *extra = nullptr) {
  if (err != cudaSuccess) {
    err = cudaGetLastError();
    fprintf(stderr, "!!! CUDA Error in operator: %s:%d %s\n", file, line, cudaGetErrorString(err));
  }
}

#define CudaCheckLaunch(f) _cudaCheckKernelLaunch((f), )
void _cudaCheckKernelLaunch(dim3 grid, dim3 block, const char *file, int line) {
  cudaError_t err = cudaPeekAtLastError();
  if (err != cudaSuccess)
    fprintf(stderr, "!!! CUDA Error after kernel launch: %s:%d %s\n", file, line, cudaGetErrorString(err));
}

#define CUDA_MAYBE_REALLOC_STATIC_TEMP(NAME, NBYTES)                           \
  if ((void *)NAME == 0 || (NBYTES) != NAME->nbytes) {                         \
    if ((void *)NAME != 0) {                                                   \
                                                                               \
      {                                                                        \
        printf("resizing NAME from %d to %d bytes\n", NAME->nbytes, (nbytes)); \
        free0(NAME);                                                           \
        NAME = nullptr;                                                        \
      }                                                                        \
    }                                                                          \
    NAME = _cudaArrayAlloc((NBYTES));                                          \
  }

#define CUDA_TMP_ALLOC(name, nbytes)                                           \
  if ((void *)name == 0) {                                                     \
    printf("allocating %d bytes for " #name "r20_vec\n", (nbytes));            \
    CudaChecked(cudaMalloc((void **)(&name), (nbytes)));                       \
  }

#include <map>

typedef struct _tempStorage {
  void *ptr;
  size_t size;
} TempStorage;

class CudaDeviceTempStorage {
public:
  CudaDeviceTempStorage() {}

  template <typename T>
  void ensureTempStorage(T **out, const std::string &name, size_t nbytes) {
    if (*out != nullptr) {
    }
  }

private:
  std::map<std::string, TempStorage> _storage;
};

uint64_t next_pow2(uint64_t x) {
  return (__builtin_popcount(x) == 1 || x == 1) ? x : 1 << (64 - __builtin_clzl(x - 1));
}

#define setupGrid(GRID, THREAD, X, Y, Z) _setupGrid(#GRID, GRID, THREAD, X, Y, Z)

// Heuristic thread block sizing based on the grid
// Not perfect, but it'll do for now
inline void _setupGrid(const char *gridName, dim3 &grid, dim3 &threadBlock, int x_size, int y_size,
                      int z_size) {
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

  grid = dim3((int)ceil((float)x_size / (float)threadBlock.x), (int)ceil((float)y_size / (float)threadBlock.y),
              (int)ceil((float)z_size / (float)threadBlock.z));
}

inline int _cudaGetCurrentDevice() {
    int device = -1;
    CudaChecked(cudaGetDevice(&device));
    return device;
}

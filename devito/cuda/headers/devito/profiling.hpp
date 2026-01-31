#ifndef _DEVITO_CUDA_PROFILING_H
#define _DEVITO_CUDA_PROFILING_H

#include <cuda_runtime.h>
#include <vector>

#include <devito/logging.hpp>

namespace devito {
namespace cuda {

class CudaSectionTimer {
public:
  CudaSectionTimer(cudaStream_t stream, double *section_ptr)
      : _stream(stream), _section_ptr(section_ptr), _start_event(nullptr),
        _end_event(nullptr) {}

  void start() {
    if (_start_event) {
      return;
    }
    cudaEventCreate(&_start_event);
    cudaEventRecord(_start_event, _stream);
  }

  void stop() {
    cudaEventCreate(&_end_event);
    cudaEventRecord(_end_event, _stream);
  }

  ~CudaSectionTimer() {
    if (_start_event) {
      cudaEventDestroy(_start_event);
      _start_event = nullptr;
    }
    if (_end_event) {
      cudaEventDestroy(_end_event);
      _end_event = nullptr;
    }
  }

  void resolve() {
    if (_start_event == nullptr) {
      return;
    }

    if (_end_event == nullptr) {
      debug("CudaSectionTimer: stop() was not called, possible codegen bug");
      return;
    }
    cudaEventSynchronize(_end_event);

    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, _start_event, _end_event);

    *_section_ptr +=
        static_cast<double>(milliseconds) / 1000.0; // Convert to seconds

    cudaEventDestroy(_start_event);
    cudaEventDestroy(_end_event);
    _start_event = nullptr;
    _end_event = nullptr;
  }

private:
  cudaStream_t _stream;
  double *_section_ptr;
  cudaEvent_t _start_event;
  cudaEvent_t _end_event;
};

class CudaSectionTimers {
public:
  CudaSectionTimers(cudaStream_t timer_stream) : _timer_stream(timer_stream) {
    // This should probably be enough - normally we target a few thousand time
    // steps at most, and large elastic operators may have 30-50 sections.
    _timers.reserve(131072);
  }

  ~CudaSectionTimers() { _timers.clear(); }

  CudaSectionTimer &startNewTimer(double *section_ptr) {
    _timers.emplace_back(_timer_stream, section_ptr);
    _timers.back().start();
    return _timers.back();
  }

  void resolveTimers() {
    debug("resolveTimers(): resolving %zu timers", _timers.size());
    for (auto &timer : _timers) {
      timer.resolve();
    }
  }

private:
  std::vector<CudaSectionTimer> _timers;
  cudaStream_t _timer_stream;
};

} // namespace cuda
} // namespace devito

#define CUDA_START_TIMER(T, S)                                                 \
  auto &_timer_##S = _cuda_section_timers.startNewTimer(&(T->S));

#define CUDA_STOP_TIMER(ST) _timer_##ST.stop();

#endif // _DEVITO_CUDA_PROFILING_H

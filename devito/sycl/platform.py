from devito.arch.archinfo import AmdDevice, Cpu64, NvidiaDevice

__all__ = ["SYCL_NVIDIA", "SYCL_AMD", "SYCL_CPU", "SYCL_SPR"]

SYCL_NVIDIA = NvidiaDevice("sycl-nv")
SYCL_AMD = AmdDevice("sycl-amd")
SYCL_CPU = Cpu64("sycl-cpu")
SYCL_SPR = Cpu64("sycl-spr")

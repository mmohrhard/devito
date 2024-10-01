from devito.ir.iet.nodes import Call, List
from devito.passes.iet.langbase import LangBB
from devito.sycl.platform import SYCL_AMD, SYCL_NVIDIA, SYCL_CPU, SYCL_SPR


__all__ = ["LangBB"]


class SyclBB(LangBB):
    mapper = {
        "name": "SYCL",
        "headers": ["sycl/sycl.hpp", "assert.h"],
        "global-decls": [],
        SYCL_NVIDIA: None,
        SYCL_AMD: None,
        SYCL_CPU: None,
        SYCL_SPR: None,
        "aligned": lambda i: "__attribute__((aligned(%d)))" % i,
        "host-alloc": lambda i, j, k: Call("posix_memalign", (i, j, k)),
        "host-free": lambda i: Call("free", (i,)),
        "init": lambda args: List(body=[]),
        "fini": lambda args: List(body=[]),
        "map-enter-to": lambda i, j: None,
        "map-enter-to-wait": lambda i, j, k: None,
        "map-enter-alloc": lambda i, j: None,
        "map-present": lambda i, j: None,
        "map-wait": lambda i: None,
        "map-update": lambda i, j: None,
        "map-update-host-if": lambda i, j, k: None,
        "map-update-host": lambda i, j: None,
        "map-update-host-async": lambda i, j, k: None,
        "map-update-host-async-if": lambda i, j, k, l: None,
        "map-update-device": lambda i, j: None,
        "map-update-device-async": lambda i, j, k: None,
        "map-update-device-async-if": lambda i, j, k, l: None,
        "map-release": lambda i, j: None,
        "map-release-if": lambda i, j, k: None,
        "map-exit-delete": lambda i, j: None,
        "map-exit-delete-if": lambda i, j, k: None,
    }

    @classmethod
    def _map_update_host(cls, f, imask=None, condition=None):
        return None  # CudaTransfer(f, imask, condition, CudaTransferDirection.D2H)

    @classmethod
    def _map_update_device(cls, f, imask=None, condition=None):
        return None  # CudaTransfer(f, imask, condition, CudaTransferDirection.H2D)

    @classmethod
    def _map_delete(cls, f, imask=None, devicerm=None):
        return None  # CudaDealloc(f, imask, devicerm)

    @classmethod
    def _map_release(cls, f, imask=None, devicerm=None):
        return None  # CudaDealloc(f, imask, devicerm)

    @classmethod
    def _map_update_host_async(cls, f, imask=None, qid=None, condition=None):
        return None

    # CudaTransfer(f, imask, condition, CudaTransferDirection.D2H, stream=qid)

    @classmethod
    def _map_update_device_async(cls, f, imask=None, qid=None, condition=None):
        return None

    # CudaTransfer(f, imask, condition, CudaTransferDirection.H2D, stream=qid)

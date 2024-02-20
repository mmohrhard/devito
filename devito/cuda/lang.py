import cgen as c

from devito.cuda.nodes import CudaTransferDirection
from devito.ir.iet.nodes import Conditional

from devito.arch import CUDA, NVIDIAX
from devito.ir import (
    Call,
    List,
    Block,
    ParallelIteration,
    Pragma,
    Definition,
    FindSymbols,
    Lambda,
)
from devito.passes.iet.parpragma import PragmaLangBB, PragmaTransfer
from devito.symbolics import Byref, VOID, INT

from devito.passes.iet.languages.openmp import OmpRegion, OmpIteration
from devito.passes.iet.languages.utils import make_clause_reduction
from devito.passes.iet.misc import is_on_device
from devito.symbolics.extended_sympy import CondNe
from devito.tools import filter_ordered
from devito.types import Symbol
from devito.cuda.nodes import (
    KernelStream,
    HostStream,
    # NcclStream,
    MemCopyStream,
    CudaChecked,
    CudaTransfer,
    CudaAlloc,
    CudaDealloc,
)
from devito.cuda.types import JitifyProgram
from devito.types.parallel import DeviceID


__all__ = ["CudaBB", "DeviceCudaIteration"]


class DeviceCudaIteration(ParallelIteration):
    @classmethod
    def _make_construct(cls, **kwargs):
        return "acc parallel loop"

    @classmethod
    def _make_clauses(cls, ncollapse=None, reduction=None, tile=None, **kwargs):
        clauses = []

        if ncollapse:
            clauses.append("collapse(%d)" % (ncollapse or 1))
        elif tile:
            clauses.append("tile(%s)" % ",".join(str(i) for i in tile))

        if reduction:
            clauses.append(make_clause_reduction(reduction))

        indexeds = FindSymbols("indexeds").visit(kwargs["nodes"])
        deviceptrs = filter_ordered(i.name for i in indexeds if i.function._mem_local)
        presents = filter_ordered(
            i.name
            for i in indexeds
            if (is_on_device(i, kwargs["gpu_fit"]) and i.name not in deviceptrs)
        )

        # The NVC 20.7 and 20.9 compilers have a bug which triggers data movement for
        # indirectly indexed arrays (e.g., a[b[i]]) unless a present clause is used
        if presents:
            clauses.append("present(%s)" % ",".join(presents))

        if deviceptrs:
            clauses.append("deviceptr(%s)" % ",".join(deviceptrs))

        return clauses

    @classmethod
    def _process_kwargs(cls, **kwargs):
        kwargs = super()._process_kwargs(**kwargs)

        kwargs.pop("gpu_fit", None)

        kwargs.pop("schedule", None)
        kwargs.pop("parallel", None)
        kwargs.pop("chunk_size", None)
        kwargs.pop("nthreads", None)
        kwargs.pop("tile", None)

        return kwargs


class CudaBB(PragmaLangBB):
    mapper = {
        # Misc
        "name": "CUDA",
        "headers": [
            "cuda.h",
            "cuda_runtime_api.h",
            "nvtx3/nvToolsExt.h",
            "stdio.h",
            "assert.h",
            "devito/devito_cuda.cuh",
            "nccl.h",
            "devito/jitify.hpp",
        ],
        "global-decls": [],
        # Platform mapping
        CUDA: None,
        NVIDIAX: None,
        # Runtime library
        "aligned": lambda i: "__attribute__((aligned(%d)))" % i,
        "init": lambda args: List(
            body=[
                Conditional(CondNe(DeviceID(), -1), Call("cudaSetDevice", (DeviceID(),))),
                Call("ENSURE_STREAM", (HostStream(),)),
                Call("ENSURE_STREAM", (MemCopyStream(),)),
                Call("ENSURE_STREAM", (KernelStream(),)),
                # Call("ENSURE_STREAM", (NcclStream(),)),
                Call("ENSURE_CACHE", ()),
                Definition(
                    JitifyProgram("program"),
                    initvalue=Call("kernel_cache.program", ("_cudaKernels", 0)),
                ),
                Call("nvtxRangePush", ("__FUNCTION__",)),
            ]
        ),
        "fini": lambda args: List(body=[Call("nvtxRangePop")]),
        "num-devices": lambda args, retobj: Block(
            body=[
                c.Initializer(c.Value("int", "_num_devices"), 0),
                Call("cudaGetDeviceCount", (INT(Byref(retobj), "*"))),
            ]
        ),
        "set-device": lambda device: CudaChecked(Call("cudaSetDevice", (device,))),
        # Pragmas
        "atomic": None,  # CUDA doesn't use a pragma for this
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
        "memcpy-to-device": lambda i, j, k: Call("acc_memcpy_to_device", [i, j, k]),
        "memcpy-to-device-wait": lambda i, j, k, l: Lambda(
            body=[
                Call("acc_memcpy_to_device_async", [i, j, k, l]),
                Call("acc_wait", [l]),
            ]
        ),
        "device-get":
        # calls a helper function since we expect a return value
        Call("_cudaGetCurrentDevice"),
        "device-alloc": lambda i, *a, retobj=None: CudaChecked(
            Call(
                "cudaMalloc",
                (
                    VOID(Byref(retobj), "**"),
                    i,
                ),
            )
        ),
        "device-free": lambda i, *a: CudaChecked(Call("cudaFree", (i,))),
        "host-alloc": lambda i, j, k:
        # this isn't really 'host', it's 'high bandwidth memory'
        CudaChecked(
            Call(
                "cudaMallocHost",
                (
                    i,
                    k,
                ),
            )
        ),
        "host-free": lambda i: CudaChecked(Call("cudaFreeHost", (i,))),
        "wait-event": lambda i, j: CudaChecked(
            Call(
                "cudaStreamWaitEvent",
                (
                    i,
                    j,
                ),
            )
        ),
        "create-event": lambda i: CudaChecked(
            Call("cudaEventCreateWithFlags", (Byref(i), "cudaEventDisableTiming"))
        ),
        "destroy-event": lambda i: CudaChecked(Call("cudaEventDestroy", (i,))),
        "record-event": lambda i, j: CudaChecked(Call("cudaEventRecord", (i, j))),
    }

    Region = OmpRegion
    HostIteration = OmpIteration  # Host parallelism still goes via OpenMP
    DeviceIteration = DeviceCudaIteration

    @classmethod
    def _map_alloc(cls, f, imask=None, condition=None):
        return CudaAlloc(f, imask, condition)

    @classmethod
    def _map_to_wait(cls, f, imask=None, qid=None):
        return None

    @classmethod
    def _map_present(cls, f, imask=None):
        return PragmaTransfer(cls.mapper["map-present"], f, imask)

    @classmethod
    def _map_wait(cls, qid=None):
        return Pragma(cls.mapper["map-wait"], qid)

    @classmethod
    def _map_update_host(cls, f, imask=None, condition=None):
        return CudaTransfer(f, imask, condition, CudaTransferDirection.D2H)

    @classmethod
    def _map_update_device(cls, f, imask=None, condition=None):
        return CudaTransfer(f, imask, condition, CudaTransferDirection.H2D)

    @classmethod
    def _map_delete(cls, f, imask=None, devicerm=None):
        return CudaDealloc(f, imask, devicerm)

    @classmethod
    def _map_release(cls, f, imask=None, devicerm=None):
        return CudaDealloc(f, imask, devicerm)

    @classmethod
    def _map_update_host_async(cls, f, imask=None, qid=None, condition=None):
        return CudaTransfer(f, imask, condition, CudaTransferDirection.D2H, stream=qid)

    @classmethod
    def _map_update_device_async(cls, f, imask=None, qid=None, condition=None):
        return CudaTransfer(f, imask, condition, CudaTransferDirection.H2D, stream=qid)

    @classmethod
    def _map_wait_event(cls, e, stream=None):
        stream = stream if stream is not None else 0
        return List(body=[cls.mapper["wait-event"](stream, e.handle)])

    @classmethod
    def _map_wait_recreate_event(cls, e, stream=None):
        stream = stream if stream is not None else 0
        return List(
            body=[
                cls.mapper["wait-event"](stream, e.handle),
                cls.mapper["destroy-event"](e.handle),
                cls.mapper["create-event"](e.handle),
            ]
        )

    @classmethod
    def _map_create_event(cls, e):
        return cls.mapper["create-event"](e.handle)

    @classmethod
    def _map_recreate_event(cls, e):
        return List(
            body=[
                cls.mapper["destroy-event"](e.handle),
                cls.mapper["create-event"](e.handle),
            ]
        )

    @classmethod
    def _map_fire_event(cls, e, stream=None):
        stream = stream if stream is not None else 0
        return cls.mapper["record-event"](e.handle, stream)

    @classmethod
    def _get_num_devices(cls, platform):
        ngpus = Symbol(name="_num_gpus")
        return ngpus, List(
            body=[
                c.Initializer(c.Value("int", "_num_gpus"), 0),
                Call("cudaGetDeviceCount", (INT(Byref(ngpus), "*"))),
            ]
        )

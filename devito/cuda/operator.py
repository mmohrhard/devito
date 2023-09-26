from functools import partial

from devito.core.gpu import (
    DeviceNoopOperator,
    DeviceAdvOperator,
    DeviceFsgOperator,
    DeviceOperatorMixin,
    make_callbacks,
)
from devito.core.operator import CustomOperator
from devito.cuda.target import DeviceCudaTarget

from devito.passes.equations import collect_derivatives
from devito.passes.clusters import (
    Lift,
    Streaming,
    Tasker,
    blocking,
    buffering,
    cire,
    cse,
    factorize,
    fission,
    fuse,
    optimize_pows,
)
from devito.passes.iet import mpiize, hoist_prodders, is_on_device

__all__ = [
    "DeviceNoopCudaOperator",
    "DeviceAdvCudaOperator",
    "DeviceAdvCudaOperator",
    "DeviceFsgCudaOperator",
    "DeviceCustomCudaOperator",
]
# CUDA


class DeviceCudaOperatorMixin(object):
    from devito.cuda.codegen import CudaCGen

    _Target = DeviceCudaTarget
    _CodeGen = CudaCGen

    @classmethod
    def _normalize_kwargs(cls, **kwargs):
        oo = kwargs["options"]
        oo.pop("openmp", None)

        kwargs = super()._normalize_kwargs(**kwargs)
        oo["cuda"] = True

        return kwargs


class DeviceNoopCudaOperator(DeviceCudaOperatorMixin, DeviceNoopOperator):
    pass


class DeviceAdvCudaOperator(DeviceCudaOperatorMixin, DeviceAdvOperator):
    pass


class DeviceFsgCudaOperator(DeviceCudaOperatorMixin, DeviceFsgOperator):
    pass


class DeviceCustomCudaOperator(
    DeviceCudaOperatorMixin, DeviceOperatorMixin, CustomOperator
):
    @classmethod
    def _make_dsl_passes_mapper(cls, **kwargs):
        return {
            "collect-derivs": collect_derivatives,
        }

    @classmethod
    def _make_iet_passes_mapper(cls, **kwargs):
        from devito.cuda.passes import cuda_eventify, cuda_linearize

        options = kwargs["options"]
        platform = kwargs["platform"]
        compiler = kwargs["compiler"]
        sregistry = kwargs["sregistry"]

        parizer = cls._Target.Parizer(sregistry, options, platform, compiler)
        orchestrator = cls._Target.Orchestrator(sregistry)

        return {
            "parallel": parizer.make_parallel,
            "cuda": parizer.make_parallel,
            "orchestrate": partial(orchestrator.process),
            "pthreadify": partial(cuda_eventify, sregistry=sregistry),
            "mpi": partial(mpiize, **kwargs),
            "linearize": partial(
                cuda_linearize, mode=options.get("linearize", None), sregistry=sregistry
            ),
            "prodders": partial(hoist_prodders),
            "init": partial(parizer.initialize, options=options),
        }

    @classmethod
    def _make_clusters_passes_mapper(cls, **kwargs):
        from devito.cuda.passes import cuda_memcpy

        options = kwargs["options"]
        platform = kwargs["platform"]
        sregistry = kwargs["sregistry"]

        # Callbacks used by `Tasking` and `Streaming`
        runs_on_host, reads_if_on_host = make_callbacks(options)

        # Callback used by `buffering`
        def callback(f):
            if not is_on_device(f, options["gpu-fit"]):
                return [f.time_dim]
            else:
                return None

        return {
            "buffering": lambda i: buffering(i, callback, sregistry, options),
            "blocking": lambda i: blocking(i, sregistry, options),
            "cuda-memcpy": lambda i: cuda_memcpy(i, sregistry=sregistry),
            "tasking": Tasker(runs_on_host, sregistry).process,
            "streaming": Streaming(reads_if_on_host, sregistry).process,
            "factorize": factorize,
            "fission": fission,
            "fuse": lambda i: fuse(i, options=options),
            "lift": lambda i: Lift().process(
                cire(i, "invariants", sregistry, options, platform)
            ),
            "cire-sops": lambda i: cire(i, "sops", sregistry, options, platform),
            "cse": lambda i: cse(i, sregistry),
            "opt-pows": optimize_pows,
            "topofuse": lambda i: fuse(i, toposort=True, options=options),
        }

    _known_passes = (
        # DSL
        "collect-derivs",
        # Expressions
        "buffering",
        # Clusters
        "blocking",
        "tasking",
        "streaming",
        "factorize",
        "fission",
        "fuse",
        "lift",
        "cire-sops",
        "cse",
        "opt-pows",
        "topofuse",
        # IET
        "orchestrate",
        "pthreadify",
        "parallel",
        "mpi",
        "linearize",
        "prodders",
        # CUDA
        "cuda",
        "cuda-events",
        "cuda-memcpy",
    )
    _known_passes_disabled = ("denormals", "simd")
    assert not (set(_known_passes) & set(_known_passes_disabled))

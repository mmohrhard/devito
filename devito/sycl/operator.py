from functools import partial

from devito.core.gpu import (
    DeviceNoopOperator,
    DeviceAdvOperator,
    DeviceFsgOperator,
    DeviceOperatorMixin,
    make_callbacks,
)
from devito.core.operator import CustomOperator, ParTile
from devito.exceptions import InvalidOperator
from devito.passes.iet.misc import relax_incr_dimensions
from devito.sycl.passes.lowering import lower_sycl
from devito.sycl.target import DeviceSyclTarget

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
from devito.tools.timing import timed_pass
import numpy as np

from devito.tools.utils import as_tuple

__all__ = [
    "DeviceNoopSyclOperator",
    "DeviceAdvSyclOperator",
    "DeviceAdvSyclOperator",
    "DeviceFsgSyclOperator",
    "DeviceCustomSyclOperator",
]
# CUDA


class DeviceSyclOperatorMixin(object):
    from devito.sycl.codegen import SyclCGen

    _Target = DeviceSyclTarget
    _CodeGen = SyclCGen

    @classmethod
    def _normalize_kwargs(cls, **kwargs):
        o = {}
        oo = kwargs["options"]

        # Execution modes
        o["mpi"] = oo.pop("mpi")
        o["parallel"] = True

        # Buffering
        o["buf-async-degree"] = oo.pop("buf-async-degree", None)

        # Fusion
        o["fuse-tasks"] = oo.pop("fuse-tasks", False)

        # Blocking
        o["blockinner"] = oo.pop("blockinner", True)
        o["blocklevels"] = oo.pop("blocklevels", cls.BLOCK_LEVELS)
        o["blockeager"] = oo.pop("blockeager", cls.BLOCK_EAGER)
        o["blocklazy"] = oo.pop("blocklazy", not o["blockeager"])
        o["blockrelax"] = oo.pop("blockrelax", cls.BLOCK_RELAX)
        o["skewing"] = oo.pop("skewing", False)

        # CIRE
        o["min-storage"] = False
        o["cire-rotate"] = False
        o["cire-maxpar"] = oo.pop("cire-maxpar", True)
        o["cire-ftemps"] = oo.pop("cire-ftemps", False)
        o["cire-mingain"] = oo.pop("cire-mingain", cls.CIRE_MINGAIN)
        o["cire-schedule"] = oo.pop("cire-schedule", cls.CIRE_SCHEDULE)
        o["cire-aggressive"] = oo.pop("cire-aggressive", False)

        # GPU parallelism
        o["par-tile"] = ParTile(oo.pop("par-tile", False), default=(32, 4))
        o["par-collapse-ncores"] = 1  # Always collapse (meaningful if `par-tile=False`)
        o["par-collapse-work"] = 1  # Always collapse (meaningful if `par-tile=False`)
        o["par-chunk-nonaffine"] = oo.pop("par-chunk-nonaffine", cls.PAR_CHUNK_NONAFFINE)
        o["par-dynamic-work"] = np.inf  # Always use static scheduling
        o["par-nested"] = np.inf  # Never use nested parallelism
        o["par-disabled"] = oo.pop("par-disabled", True)  # No host parallelism by default
        o["gpu-fit"] = as_tuple(oo.pop("gpu-fit", cls._normalize_gpu_fit(**kwargs)))
        o["gpu-nofit"] = as_tuple(oo.pop("gpu-nofit", None))

        # Misc
        o["optcomms"] = oo.pop("optcomms", True)
        oo.pop("linearize", None)
        oo.pop("openmp", None)
        o["linearize"] = True
        o["lower_sycl"] = oo.pop("lower_sycl", {})
        o["mapify-reduce"] = oo.pop("mapify-reduce", cls.MAPIFY_REDUCE)
        o["realign_iet"] = oo.pop("realign_iet", True)

        if oo:
            raise InvalidOperator(
                "Unsupported optimization options: [%s]" % ", ".join(list(oo))
            )

        kwargs["options"].update(o)

        return kwargs

    @classmethod
    def _make_iet_passes_mapper(cls, **kwargs):
        from devito.sycl.passes import sycl_linearize, lower_sycl

        options = kwargs["options"]
        platform = kwargs["platform"]
        compiler = kwargs["compiler"]
        sregistry = kwargs["sregistry"]

        parizer = cls._Target.Parizer(sregistry, options, platform, compiler)
        orchestrator = cls._Target.Orchestrator(sregistry)

        return {
            "parallel": parizer.make_parallel,
            "sycl": parizer.make_parallel,
            "orchestrate": partial(orchestrator.process),
            # "pthreadify": partial(sycl_eventify, sregistry=sregistry),
            "mpi": partial(mpiize, **kwargs),
            "linearize": partial(sycl_linearize, mode=True, sregistry=sregistry),
            "prodders": partial(hoist_prodders),
            "init": partial(parizer.initialize, options=options),
            "lower_sycl": partial(lower_sycl),
        }


class DeviceNoopSyclOperator(DeviceSyclOperatorMixin, DeviceNoopOperator):
    pass


class DeviceAdvSyclOperator(DeviceSyclOperatorMixin, DeviceAdvOperator):
    @classmethod
    @timed_pass(name="specializing.IET")
    def _specialize_iet(cls, graph, **kwargs):
        options = kwargs["options"]
        platform = kwargs["platform"]
        compiler = kwargs["compiler"]
        sregistry = kwargs["sregistry"]

        # Distributed-memory parallelism
        mpiize(graph, **kwargs)

        # Lower BlockDimensions so that blocks of arbitrary shape may be used
        relax_incr_dimensions(graph)

        # GPU parallelism
        parizer = cls._Target.Parizer(sregistry, options, platform, compiler)
        parizer.make_parallel(graph)
        parizer.initialize(graph, options=options)

        # Misc optimizations
        hoist_prodders(graph)

        # Symbol definitions
        cls._Target.DataManager(sregistry, options).process(graph)

        # Linearize n-dimensional Indexeds
        linearizer = cls._Target.Linearizer
        linearizer(graph, mode=options["linearize"], sregistry=sregistry)

        lower_sycl(graph, mode=options["lower_sycl"], sregistry=sregistry)
        return graph


class DeviceFsgSyclOperator(DeviceSyclOperatorMixin, DeviceFsgOperator):
    pass


class DeviceCustomSyclOperator(
    DeviceSyclOperatorMixin, DeviceOperatorMixin, CustomOperator
):
    @classmethod
    def _make_dsl_passes_mapper(cls, **kwargs):
        return {
            "collect-derivs": collect_derivatives,
        }

    @classmethod
    def _make_clusters_passes_mapper(cls, **kwargs):
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
        # SYCL
        "sycl",
        "lower_sycl",
    )
    _known_passes_disabled = ("denormals", "simd")
    assert not (set(_known_passes) & set(_known_passes_disabled))

from devito.arch.archinfo import platform_registry
from .target import *  # noqa
from .operator import *  # noqa
from .platform import *  # noqa

from devito.operator.registry import operator_registry


operator_registry._languages = ("C", "openmp", "openacc", "cuda", "sycl")

platform_registry["sycl-nv"] = SYCL_NVIDIA
platform_registry["sycl-amd"] = SYCL_AMD
platform_registry["sycl-cpu"] = SYCL_CPU
platform_registry["sycl-spr"] = SYCL_SPR
for platform in (SYCL_NVIDIA, SYCL_AMD, SYCL_CPU, SYCL_SPR):
    operator_registry.add(DeviceCustomSyclOperator, platform.__class__, "custom", "sycl")
    operator_registry.add(DeviceNoopSyclOperator, platform.__class__, "noop", "sycl")
    operator_registry.add(DeviceAdvSyclOperator, platform.__class__, "advanced", "sycl")
    operator_registry.add(
        DeviceFsgSyclOperator, platform.__class__, "advanced-fsg", "sycl"
    )

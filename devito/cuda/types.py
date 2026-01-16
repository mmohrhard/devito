import ctypes

import sympy

from devito.symbolics.extended_sympy import ReservedWord
from devito.types import Global, Scalar
from devito.types.array import ArrayMapped
from devito.types.parallel import DeviceSymbol
from devito.types.utils import DimensionTuple

__all__ = [
    "CudaError",
    "CudaEvent",
    "CudaStream",
    "NullPointer",
    "JitifyCache",
    "JitifyProgram",
    "EnlargedBuffer",
]


class cudaEvent_t(ctypes.Structure):
    pass


class cudaStream_t(ctypes.Structure):
    pass


c_cudaEvent_p = ctypes.POINTER(cudaEvent_t)


class CudaError(Scalar):
    def __init__(self, name):
        super().__init__(name=name, dtype=ctypes.c_int)

    @property
    def _C_typename(self):
        return "cudaError_t"


class CudaEvent(Scalar):
    def __init__(self, name):
        super().__init__(name=name, dtype=ctypes.c_void_p)

    @property
    def _C_typename(self):
        return "cudaEvent_t"


class CudaStream(Scalar):
    def __init__(self, name):
        super().__init__(name=name, dtype=ctypes.c_void_p)

    @property
    def _C_typename(self):
        return "cudaStream_t"


class NullPointer(ReservedWord):
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls, "nullptr")


class JitifyCache(Global):
    @property
    def _C_typename(self):
        return "jitify::JitCache"

    def __init__(self, name, *args, **kwargs):
        super().__init__(name, dtype=ctypes.c_void_p)

    def __new__(cls, *args, **kwargs):
        return super().__new__(cls, "kernel_cache")


class JitifyProgram(Global):
    @property
    def _C_typename(self):
        return "jitify::Program"

    def __init__(cls, name, *args, **kwargs):
        super().__init__(name, dtype=ctypes.c_void_p)

    def __new__(cls, name, *args, **kwargs):
        return super().__new__(cls, name)


class JitifyKernelInstantiation(Scalar):
    def __init__(self, name, *_):
        super().__init__(name=name, dtype=ctypes.c_void_p)

    @property
    def _C_typename(self):
        return "jitify::KernelInstantiation"

    @property
    def _C_arg_typename(self):
        return "jitify::KernelInstantiation&"


class JitKernelTuningParams(Scalar):
    def __init__(self, name, *_):
        super().__init__(name=name, dtype=ctypes.c_void_p)

    @property
    def _C_typename(self):
        return "std::tuple<dim3, dim3> "


class EnlargedBuffer(ArrayMapped):
    def __init_finalize__(self, **kwargs):
        super().__init_finalize__(**kwargs)
        self.adjusted_dimensions = [
            d.parent if d.is_Sub else d for d in kwargs["dimensions"]
        ]

    @property
    def symbolic_shape(self):
        """
        The symbolic shape of the object. This includes the domain, halo, and
        padding regions. While halo and padding are known quantities (integers),
        the domain size is given as a symbol.
        """
        halo = [sympy.Add(*i, evaluate=False) for i in self._size_halo]
        padding = [sympy.Add(*i, evaluate=False) for i in self._size_padding]
        domain = [i.symbolic_size for i in self.adjusted_dimensions]
        ret = tuple(sympy.Add(i, j, k) for i, j, k in zip(domain, halo, padding))
        return DimensionTuple(*ret, getters=self.dimensions)


class JitOnly(DeviceSymbol):
    name = "jitonly"

    @property
    def default_value(self):
        return 0

    def _arg_values(self, **kwargs):
        try:
            # Enforce 1 or 0
            return {self.name: int(bool(kwargs[self.name]))}
        except KeyError:
            return self._arg_defaults()

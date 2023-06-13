import ctypes
from devito.types import Scalar, Global
from devito.symbolics.extended_sympy import ReservedWord

__all__ = ['CudaEvent', 'CudaStream', 'NullPointer', 'JitifyCache', 'JitifyProgram']

class cudaEvent_t(ctypes.Structure):
    pass

class cudaStream_t(ctypes.Structure):
    pass

c_cudaEvent_p = ctypes.POINTER(cudaEvent_t)

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
    def __init__(cls, name, *args, **kwargs):
        super().__init__(name, dtype=ctypes.c_void_p)

    def __new__(cls, *args):
        return super().__new__(cls, "kernel_cache")

    def __new__(cls, *args, **kwargs):
        return super().__new__(cls, "kernel_cache")
    
class JitifyProgram(Global):
    @property
    def _C_typename(self):
        return "jitify::Program"    
    def __init__(cls, name, *args, **kwargs):
        super().__init__(name, dtype=ctypes.c_void_p)

    def __new__(cls, name, *args):
        return super().__new__(cls, name)

    def __new__(cls, name, *args, **kwargs):
        return super().__new__(cls, name)
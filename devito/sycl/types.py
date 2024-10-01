import ctypes


from devito.ir.equations.equation import DummyEq
from devito.ir.iet.nodes import Expression
from devito.symbolics.extended_sympy import MacroArgument
from devito.tools.utils import ctypes_to_cstr, dtype_to_ctype
from devito.types.basic import Scalar
from devito.types.misc import Global

__all__ = [
    "SyclQueue",
    "SyclEvent",
    "SyclHandler",
    "SyclSpecializationConstantWrite",
    "SyclSpecializationConstant",
]


class SyclQueue(Scalar):
    def __init__(self, *args, **kwargs):
        super().__init__(self, name="sycl_queue", dtype=ctypes.c_void_p, is_const=True)

    def __new__(cls, *args, **kwargs):
        return super().__new__(
            cls, name="sycl_queue", dtype=ctypes.c_void_p, is_const=True
        )

    @property
    def _C_typename(self):
        return "sycl::queue *"


class SyclEvent(Scalar):
    def __init__(self, name):
        super().__init__(name=name, dtype=ctypes.c_void_p)

    @property
    def _C_typename(self):
        return "sycl::event"


class SyclHandler(Scalar):
    def __init__(self, name):
        super().__init__(name=name, dtype=ctypes.c_void_p)

    @property
    def _C_typename(self):
        return "sycl::handler"


class SyclSpecializationConstantWrite(Expression):
    def __init__(self, lhs, rhs):
        super().__init__(DummyEq(lhs, rhs), init=True)


class SyclSpecializationConstant(MacroArgument):
    @property
    def _C_typename(self):
        return "sycl::specialization_id<%s>" % ctypes_to_cstr(dtype_to_ctype(self.dtype))

    def __str__(self):
        return "_handler.get_specialization_constant<%s>()" % self.name

    __repr__ = __str__


class SyclSpecializationConstantDefinition(Global):
    def __init__(self, name, dtype):
        super().__init__(self, name=name, dtype=dtype, is_const=True)

    @property
    def _C_typename(self):
        return "sycl::specialization_id<%s>" % ctypes_to_cstr(dtype_to_ctype(self.dtype))

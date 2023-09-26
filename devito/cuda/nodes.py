from enum import Enum
from cached_property import cached_property

from devito.data import FULL
from devito.tools.utils import as_tuple, filter_ordered, flatten

from devito.ir.iet.nodes import (
    Call,
    CallableBody,
    DeviceCall,
    DeviceFunction,
    Expression,
    ExprStmt,
    Node,
    Global,
    Definition,
    Transfer,
    CLiteral,
)

from devito.ir.iet.efunc import AsyncCall, AsyncCallable
from devito.symbolics import ccode
from devito.types import Scalar
from devito.ir.equations import DummyEq
from devito.cuda.types import CudaStream

import ctypes as c

__all__ = [
    "CudaCall",
    "CudaCallable",
    "CudaCallableBody",
    "CudaTransferDirection",
    "CudaConstantWrite",
    "CudaConstantDecl",
    "CudaKernelPointerCast",
    "TemplateParameter",
    "CudaKernelTuner",
    "CudaTunedKernel",
]


class CudaCall(DeviceCall):
    is_Call = True

    def __init__(
        self,
        name=None,
        grid=None,
        threads=None,
        preferred_block=None,
        preferred_sub_block=None,
        template_arguments=None,
        arguments=None,
        kernel=None,
        writes=None,
        types=None,
        stream=None,
    ):
        super().__init__(name, arguments, None, writes=writes, types=types)
        self._grid = grid
        self._threads = threads
        self._stream = stream
        self._template_arguments = template_arguments
        self._kernel = kernel
        self._preferred_block = preferred_block
        self._preferred_sub_block = preferred_sub_block

    @property
    def grid(self):
        return self._grid

    @property
    def threads(self):
        return self._threads

    @property
    def stream(self):
        return self._stream

    @property
    def preferred_block(self):
        return self._preferred_block

    @property
    def preferred_sub_block(self):
        return self._preferred_sub_block

    @property
    def template_arguments(self):
        return self._template_arguments

    def __repr__(self):
        ret = "" if self.retobj is None else "%s = " % self.retobj
        return "%sCudaCall::%s(...)" % (ret, self.name)


class CudaCallableBody(CallableBody):
    _traversable = [
        "unpacks",
        "casts",
        "init",
        "allocs",
        "maps",
        "objs",
        "body",
        "unmaps",
        "frees",
    ]

    def __init__(
        self,
        body,
        init=None,
        unpacks=None,
        allocs=None,
        casts=None,
        objs=None,
        maps=None,
        unmaps=None,
        frees=None,
        fini=None,
    ):
        super().__init__(
            body, init, unpacks, allocs, casts, objs, maps, unmaps, frees, fini
        )

    def __repr__(self):
        return (
            "<CudaCallableBody <unpacks=%d, allocs=%d, casts=%d, maps=%d, "
            "objs=%d> <unmaps=%d, frees=%d>>"
            % (
                len(self.unpacks),
                len(self.allocs),
                len(self.casts),
                len(self.maps),
                len(self.objs),
                len(self.unmaps),
                len(self.frees),
            )
        )

    @property
    def used_globals(self):
        from devito.ir.iet.visitors import FindNodes

        return filter_ordered(FindNodes(Global).visit(self.body))


class CudaCallable(DeviceFunction):
    is_Callable = True

    _traversable = ["body"]

    _defines = None

    def __init__(
        self,
        name=None,
        body=None,
        parameters=None,
        defines=None,
        template_parameters=None,
        preferred_block=None,
        preferred_sub_block=None,
        block_dims=None,
    ):
        if isinstance(body, CallableBody):
            super().__init__(
                name,
                CudaCallableBody(
                    body.body,
                    body.init,
                    body.unpacks,
                    body.allocs,
                    body.casts,
                    body.objs,
                    body.maps,
                    body.unmaps,
                    body.frees,
                    body.fini,
                ),
                "void",
                parameters=parameters,
                prefix="__global__",
            )
        else:
            super().__init__(
                name,
                CudaCallableBody(body),
                "void",
                parameters=parameters,
                prefix="__global__",
            )
        self._defines = defines
        self._template_parameters = template_parameters or []
        self._preferred_block = preferred_block or (64,)
        self._preferred_sub_block = preferred_sub_block or as_tuple(
            [1] * len(self._preferred_block)
        )

    @property
    def template_parameters(self):
        return filter_ordered(
            flatten(
                self.block_parameters
                + self.sub_block_parameters
                + [x.free_symbols for x in self.body.casts]
                + self._template_parameters
            )
        )

    @property
    def preferred_block(self):
        return self._preferred_block

    @property
    def preferred_sub_block(self):
        return self._preferred_sub_block

    @property
    def block_parameters(self):
        return [x.lhs for x in self.block_arguments]

    @property
    def sub_block_parameters(self):
        return [x.lhs for x in self.sub_block_arguments]

    @property
    def block_arguments(self):
        dims = ["x", "y", "z"]
        return [
            DummyEq(TemplateParameter(f"_block_{d[1]}"), d[0])
            for d in zip(self._preferred_block, dims)
        ]

    @property
    def sub_block_arguments(self):
        dims = ["x", "y", "z"]
        return [
            DummyEq(TemplateParameter(f"_sub_block_{d[1]}"), d[0])
            for d in zip(self._preferred_sub_block, dims)
        ]

    @property
    def template_arguments(self):
        return filter_ordered(flatten([x.expr_symbols for x in self.body.casts]))

    @property
    def defines(self):
        return [
            x
            for x in flatten([self.parameters, self._defines, self._template_parameters])
            if x
        ]

    @cached_property
    def writes(self):
        from devito.ir.iet.visitors import FindNodes

        return filter_ordered(
            flatten([x.write for x in FindNodes(Expression).visit(self)])
        )


class CudaConstantWrite(Expression):
    def __init__(self, symbol, rhs, **kwargs):
        super().__init__(DummyEq(symbol, rhs))
        self.lhs = symbol
        self.rhs = rhs

    @property
    def free_symbols(self):
        return []

    def writes(self):
        return self.lhs

    def __repr__(self):
        return "<CudaConstant (%s)=%s>" % (self.lhs, self.rhs)


class CudaConstantDecl(Definition):
    def __init__(self, function, dtype):
        super().__init__(function)
        self.dtype = dtype

    def __repr__(self):
        return "<CudaConstantDecl(%s)>" % self.function


class CudaTransferDirection(Enum):
    H2D = 1
    D2H = 2


class TemplateParameter(Global):
    pass


class CudaKernelPointerCast(ExprStmt, Node):
    """
    A node encapsulating a cast of a raw void pointer to a non-void,
    potentially multi-dimensional array.
    """

    is_PointerCast = True

    def __init__(self, function, obj=None, alignment=True, flat=None):
        self.function = function
        self.obj = obj
        self.alignment = alignment
        self.flat = flat

    def __repr__(self):
        return "<CudaKernelPointerCast(%s)>" % self.function

    @property
    def castshape(self):
        """
        The shape used in the left-hand side and right-hand side
        of the CudaKernelPointerCast.
        """
        if self.function.is_ArrayBasic:
            calc = lambda x: self.function.symbolic_shape[x]
        else:
            calc = lambda x: self.function._C_get_field(FULL, x).size

        return tuple(
            DummyEq(TemplateParameter("%s_sz_%s" % (self.function.name, d.name)), calc(d))
            for d in self.function.dimensions[1:]
        )

    @property
    def functions(self):
        return (self.function,)

    @property
    def expr_symbols(self):
        return self.castshape

    @property
    def defines(self):
        return (self.function.indexed,)

    @property
    def free_symbols(self):
        return [x.lhs for x in self.castshape]


class CudaTunedKernel(Scalar):
    """
    A symbol containing a tuned CUDA kernel ready for execution.
    """

    def __init__(cls, name):
        super().__init__(name=name, dtype=c.c_void_p)

    def _C_typename(self):
        return "KernelInstantiation&"


class CudaKernelTuner(ExprStmt, Node):
    """
    A node encapsulating a just-in-time tuning operation for a CUDA kernel given
    the known input data shapes.
    """

    def __init__(self, output_kernel=None):
        self._output_kernel = output_kernel

    def __repr__(self):
        return "<CudaKernelTuner(%s)>" % "a"

    @property
    def defines(self):
        return (self._output_kernel,)


class KernelStream(CudaStream, Global):
    def __init__(cls, *args, **kwargs):
        super().__init__("kernel_stream")

    def __new__(cls, *args, **kwargs):
        return super().__new__(cls, "kernel_stream")


class MemCopyStream(CudaStream, Global):
    def __init__(cls, *args, **kwargs):
        super().__init__("memcpy_stream")

    def __new__(cls, *args, **kwargs):
        return super().__new__(cls, "memcpy_stream")


class HostStream(CudaStream, Global):
    def __init__(cls, *args, **kwargs):
        super().__init__("host_stream")

    def __new__(cls, *args, **kwargs):
        return super().__new__(cls, "host_stream")


class NcclStream(CudaStream, Global):
    def __init__(cls, *args, **kwargs):
        super().__init__("nccl_stream")

    def __new__(cls, *args, **kwargs):
        return super().__new__(cls, "nccl_stream")


class CudaStorage:
    def __init__(self, function, imask=None):
        self._imask = imask
        self._function = function

    @property
    def function(self):
        return self._function

    @cached_property
    def device_storage(self):
        return "%s->%s" % (self.function._C_name, self.function._C_field_device_data)

    @cached_property
    def host_storage(self):
        return "%s->%s" % (self.function._C_name, self.function._C_field_data)

    @cached_property
    def operator_allocated(self):
        return "%s->%s" % (
            self.function._C_name,
            self.function._C_field_operator_allocated,
        )

    @cached_property
    def size(self):
        return ("sizeof(%s) * " % (self.function.indexed._C_typedata)) + " * ".join(
            "(" + ccode(j) + ")" for i, j in self.sections
        )

    @property
    def imask(self):
        return self._imask

    @cached_property
    def name(self):
        return self.function._C_name

    @cached_property
    def sections(self):
        from devito.passes.iet.langbase import make_sections_from_imask

        return make_sections_from_imask(self.function, self.imask)


class CudaTransfer(CudaStorage, Transfer, Node):
    """
    A data transfer between host and CUDA device.
    """

    def __init__(
        self,
        function,
        imask=None,
        condition=None,
        direction=CudaTransferDirection.H2D,
        delete=None,
        stream=None,
    ):
        super().__init__(function, imask)

        self._direction = direction
        self._condition = condition
        self._delete = delete
        self._stream = stream

    @cached_property
    def direction(self):
        return self._direction

    @property
    def condition(self):
        return self._condition

    @property
    def functions(self):
        return (self.function,)

    @property
    def delete(self):
        return self._delete

    @property
    def stream(self):
        return self._stream

    @cached_property
    def expr_symbols(self):
        retval = [self.function.indexed]
        for i in (self.condition,) + tuple(flatten(self.sections)):
            try:
                retval.extend(i.free_symbols)
            except AttributeError:
                pass
        return tuple(retval)


class CudaAlloc(CudaStorage, Node):
    def __init__(self, function, imask=None, condition=None):
        super().__init__(function, imask)
        self._condition = condition

    @property
    def condition(self):
        return self._condition

    @cached_property
    def expr_symbols(self):
        retval = [self.function.indexed]
        for i in (self.condition,) + tuple(flatten(self.sections)):
            try:
                retval.extend(i.free_symbols)
            except AttributeError:
                pass
        return tuple(retval)


class CudaDealloc(CudaStorage, Node):
    def __init__(self, function, imask=None, condition=None):
        super().__init__(function, imask)
        self._condition = condition

    @property
    def condition(self):
        return self._condition

    @cached_property
    def expr_symbols(self):
        retval = [self.function.indexed]
        for i in (self.condition,) + tuple(flatten(self.sections)):
            try:
                retval.extend(i.free_symbols)
            except AttributeError:
                pass
        return tuple(retval)


class CudaCheckError(CLiteral):
    def __init__(self):
        super().__init__(
            "if (cudaPeekAtLastError() != 0 ) { cudaError_t err "
            '= cudaGetLastError(); printf("\\n!E %s: %s\\n",'
            "cudaGetErrorName(err), cudaGetErrorString(err));}"
        )


class CudaChecked(Call):
    def __init__(self, arguments=None, name=None):
        super().__init__("CudaChecked", arguments=tuple(flatten([arguments])))

    @cached_property
    def expr_symbols(self):
        return flatten([x.expr_symbols for x in flatten(self.arguments)])


class CudaAtomicExpression(Expression):
    def __init__(self, expr, pragmas=None, init=None, operation=None):
        super().__init__(expr, pragmas, init, operation, True)


class CudaHostFuncCall(AsyncCall):
    def __init__(
        self,
        name,
        arguments=None,
        retobj=None,
        is_indirect=False,
        cast=False,
        writes=None,
        types=None,
        declares=True,
        stream=None,
    ):
        super().__init__(
            name,
            arguments=arguments,
            retobj=retobj,
            is_indirect=is_indirect,
            cast=cast,
            writes=writes,
            types=types,
            declares=declares,
        )
        self.stream = stream


class CudaHostFuncCallable(AsyncCallable):
    pass

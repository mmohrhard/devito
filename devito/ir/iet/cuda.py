from enum import Enum
from cached_property import cached_property

from devito.data import FULL
from devito.tools.utils import as_tuple, filter_ordered, flatten

from devito.ir.iet.nodes import Callable, Call, CallableBody, DummyExpr, Expression, ExprStmt, Node, Global, Definition
from devito.ir.equations import DummyEq

__all__ = ['CudaCall', 'CudaCallable', 'CudaCallableBody', 'DeviceFunction', 'DeviceCall', 'CudaTransferDirection',
           'CudaConstantWrite', 'CudaConstantDecl', 'CudaKernelPointerCast', 'TemplateParameter']

class DeviceFunction(Callable):

    """
    A Callable executed asynchronously on a device.
    """

    def __init__(self, name, body, retval='void', parameters=None, prefix='__global__'):
        super().__init__(name, body, retval, parameters=parameters, prefix=prefix)


class DeviceCall(Call):

    """
    A call to an external function executed asynchronously on a device.
    """

    pass

class CudaCall(DeviceCall):
    is_Call = True

    def __init__(self, name=None, grid=None, threads=None, preferred_block=None, preferred_sub_block=None, template_arguments=None, arguments=None, kernel=None, writes=None, types=None, stream=None):
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
    _traversable = ['unpacks', 'casts', 'init', 'allocs', 'maps', 'objs',
                    'body', 'unmaps', 'frees']
    def __init__(self, body, init=None, unpacks=None, allocs=None, casts=None,
                 objs=None, maps=None, unmaps=None, frees=None, fini=None):
        super().__init__(body, init, unpacks, allocs, casts, objs, maps, unmaps, frees, fini)

    def __repr__(self):
        return ("<CudaCallableBody <unpacks=%d, allocs=%d, casts=%d, maps=%d, "
                "objs=%d> <unmaps=%d, frees=%d>>" %
                (len(self.unpacks), len(self.allocs), len(self.casts),
                 len(self.maps), len(self.objs), len(self.unmaps),
                 len(self.frees)))

    @property
    def used_globals(self):
        from devito.ir.iet.visitors import FindNodes
        return filter_ordered(FindNodes(Global).visit(self.body))
    

class CudaCallable(DeviceFunction):
    is_Callable = True

    _traversable = ['body']

    _defines = None

    def __init__(self, name=None, body=None, parameters=None, defines=None, template_parameters=None, preferred_block=None, preferred_sub_block=None, block_dims=None):
        if isinstance(body, CallableBody):
            super().__init__(name, CudaCallableBody(body.body, body.init, body.unpacks, body.allocs, body.casts, body.objs, body.maps, body.unmaps, body.frees, body.fini), 'void', parameters=parameters, prefix='__global__')
        else:
            super().__init__(name, CudaCallableBody(body), 'void', parameters=parameters, prefix='__global__')
        self._defines = defines
        self._template_parameters = template_parameters or []
        self._preferred_block = preferred_block or (64,)
        self._preferred_sub_block = preferred_sub_block or as_tuple([1] * len(self._preferred_block))

    @property
    def template_parameters(self):
        from devito.ir.iet.visitors import FindNodes
        
        return filter_ordered(flatten(self.block_parameters + self.sub_block_parameters + [x.free_symbols for x in self.body.casts] + self._template_parameters))
    
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
        dims = ['x', 'y', 'z']
        return [DummyEq(TemplateParameter(f"_block_{d[1]}"), d[0]) for d in zip(self._preferred_block, dims)]
    
    @property
    def sub_block_arguments(self):
        dims = ['x', 'y', 'z']
        return [DummyEq(TemplateParameter(f"_sub_block_{d[1]}"), d[0]) for d in zip(self._preferred_sub_block, dims)]
    
    @property
    def template_arguments(self):
        return filter_ordered(flatten([x.expr_symbols for x in self.body.casts]))
    
    @property
    def defines(self):
        return [x for x in flatten([self.parameters, self._defines, self._template_parameters]) if x]
    
    @cached_property
    def writes(self):
        from devito.ir.iet.visitors import FindNodes
        return filter_ordered(flatten([x.write for x in FindNodes(Expression).visit(self)]))

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
        The shape used in the left-hand side and right-hand side of the CudaKernelPointerCast.
        """
        if self.function.is_ArrayBasic:
            calc = lambda x: self.function.symbolic_shape[x]
        else:
            calc = lambda x: self.function._C_get_field(FULL, x).size

        return tuple(DummyEq(TemplateParameter('%s_sz_%s' % (self.function.name, d.name)), calc(d)) for d in self.function.dimensions[1:])
        

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
     
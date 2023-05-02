from enum import Enum
from cached_property import cached_property

from devito.tools.utils import filter_ordered, flatten

from devito.ir.iet.nodes import Callable, Call, CallableBody, Expression


__all__ = ['CudaCall', 'CudaCallable', 'DeviceFunction', 'DeviceCall', 'CudaTransferDirection']

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

    def __init__(self, name=None, grid=None, threads=None, arguments=None, writes=None, types=None, stream=None):
        super().__init__(name, arguments, None, writes=writes, types=types)
        self._grid = grid
        self._threads = threads
        self._stream = stream

    @property
    def grid(self):
        return self._grid

    @property
    def threads(self):
        return self._threads
    
    @property
    def stream(self):
        return self._stream
    
    def __repr__(self):
        ret = "" if self.retobj is None else "%s = " % self.retobj
        return "%sCudaCall::%s(...)" % (ret, self.name)


class CudaCallableBody(CallableBody):
    _traversable = ['unpacks', 'init', 'allocs', 'maps', 'objs',
                    'body', 'unmaps', 'frees']
    def __init__(self, body, init=None, unpacks=None, allocs=None, casts=None,
                 objs=None, maps=None, unmaps=None, frees=None, fini=None):
        super().__init__(body, init, unpacks, allocs, None, objs, maps, unmaps, frees, fini)

    def __repr__(self):
        return ("<CudaCallableBody <unpacks=%d, allocs=%d, casts=%d, maps=%d, "
                "objs=%d> <unmaps=%d, frees=%d>>" %
                (len(self.unpacks), len(self.allocs), len(self.casts),
                 len(self.maps), len(self.objs), len(self.unmaps),
                 len(self.frees)))
    

class CudaCallable(DeviceFunction):
    is_Callable = True

    _traversable = ['body']

    _defines = None

    def __init__(self, name=None, body=None, parameters=None, defines=None):
        if isinstance(body, CallableBody):
            super().__init__(name, CudaCallableBody(body.body, body.init, body.unpacks, body.allocs, body.casts, body.objs, body.maps, body.unmaps, body.frees, body.fini), 'void', parameters=parameters, prefix='__global__')
        else:
            super().__init__(name, CudaCallableBody(body), 'void', parameters=parameters, prefix='__global__')
        self._defines = defines

    @property
    def defines(self):
        return flatten([self.parameters, self._defines])
    
    @cached_property
    def writes(self):
        from devito.ir.iet.visitors import FindNodes
        return filter_ordered(flatten([x.write for x in FindNodes(Expression).visit(self)]))

class CudaTransferDirection(Enum):
    H2D = 1
    D2H = 2
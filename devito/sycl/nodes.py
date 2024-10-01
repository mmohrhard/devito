from functools import cached_property
from devito.ir.iet.nodes import CallableBody, DeviceFunction, Expression, Node
from devito.sycl.types import SyclQueue
from devito.tools.utils import as_tuple, filter_ordered, flatten
from devito.types.misc import Global

__all__ = ["SyclKernel", "SyclKernelBody"]


class SyclKernel(DeviceFunction):
    is_Callable = True

    _traversable = ["body"]

    _defines = None

    def __init__(
        self,
        name=None,
        body=None,
        parameters=None,
        defines=None,
    ):
        if isinstance(body, CallableBody):
            super().__init__(
                name,
                SyclKernelBody(
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
            )
        else:
            super().__init__(
                name,
                SyclKernelBody(body),
                "void",
                parameters=parameters,
            )
        self._defines = defines or []

    @cached_property
    def writes(self):
        from devito.ir.iet.visitors import FindNodes

        return filter_ordered(
            flatten([x.write for x in FindNodes(Expression).visit(self)])
        )

    @property
    def defines(self):
        return as_tuple(flatten([self._defines, self.body.defines]))


class SyclKernelBody(CallableBody):
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
            "<SyclKernelBodyBody <unpacks=%d, allocs=%d, casts=%d, maps=%d, "
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


class SyclKernelLaunch(Node):
    _traversable = ["kernel", "arguments"]

    def __init__(
        self,
        name=None,
        kernel=None,
        constant_writes=None,
        dims=None,
        depends_on=None,
        queue=None,
        synchronous=True,
        **kwargs
    ):
        super().__init__(**kwargs)
        self.name = name
        self.kernel = kernel
        self.dims = dims
        self.constant_writes = constant_writes
        self.depends_on = depends_on
        self.queue = queue or SyclQueue()
        self.synchronous = synchronous

    @property
    def children(self):
        return [self.kernel]

    @property
    def expr_symbols(self):
        return as_tuple(flatten([self.kernel.expr_symbols, self.dims, self.queue]))

    def __repr__(self):
        return "<SyclKernel <%s>>" % (self.name)


class SyclAtomicExpression(Expression):
    def __init__(self, expr, pragmas=None, init=None, operation=None):
        super().__init__(expr, pragmas, init, operation, True)


class SyclAlloc(Node):
    def __init__(self, symbol=None, nbytes=0, queue=None, **kwargs):
        self.symbol = symbol
        self.nbytes = nbytes
        self.queue = queue

    @property
    def expr_symbols(self):
        return (self.symbol,)

    def __repr__(self):
        return "<SyclAlloc <%s, size=%s>>" % (self.symbol, self.nbytes)


class SyclDealloc(Node):
    def __init__(self, symbol=None, queue=None, **kwargs):
        self.symbol = symbol
        self.queue = queue

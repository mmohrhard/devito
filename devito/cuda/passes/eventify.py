import cgen as c
from collections import OrderedDict
import numpy as np

from devito.ir.iet.efunc import AsyncCall, ThreadCallable
from devito.ir.iet.nodes import (BlankLine, Call, Callable, Definition, Dereference,
                                 DummyExpr, List, PointerCast, Return)
from devito.ir.iet.visitors import FindNodes, FindSymbols, Transformer
from devito.logger import debug
from devito.passes.iet.engine import iet_pass
from devito.symbolics.extended_sympy import (VOID, Byref, FieldFromComposite,
                                             FieldFromPointer, Null, SizeOf)
from devito.symbolics.printer import ccode
from devito.tools.data_structures import Bunch, DefaultOrderedDict
from devito.tools.utils import as_list, flatten, split
from devito.types.misc import Pointer
from devito.types.parallel import QueueID, ThreadArray

from devito.cuda.nodes import CudaHostFuncCall, CudaHostFuncCallable
from devito.cuda.types import NullPointer

__all__ = ['cuda_eventify']


def cuda_eventify(graph, **kwargs):
    """
    Rewrites AsyncCalls into CUDA host launches
    """
    track = DefaultOrderedDict(lambda: Bunch(threads=None, sdata=None, extra_args=None))

    debug("performing CUDA eventification")
    lower_async_callables(graph, track=track, root=graph.root, **kwargs)
    lower_async_calls(graph, track=track, **kwargs)


class CudaSharedData(ThreadArray):

    """
    An Array of structs, each struct containing data shared by one producer and
    one consumer thread.
    """

    __rkwargs__ = list(ThreadArray.__rkwargs__) + ['cfields', 'ncfields']
    __rkwargs__.remove('fields')

    def __init_finalize__(self, *args, **kwargs):
        self.cfields = tuple(kwargs.pop('cfields', ()))
        self.ncfields = tuple(kwargs.pop('ncfields', ()))

        kwargs['fields'] = self.cfields + self.ncfields

        super().__init_finalize__(*args, **kwargs)

    @property
    def _mem_stack(self):
        return False

    @property
    def _C_name(self):
        return self.name

    @classmethod
    def __pfields_setup__(cls, **kwargs):
        fields = as_list(kwargs.get('cfields'))
        fields.extend(as_list(kwargs.get('ncfields')))
        return [(i._C_name, i._C_ctype) for i in fields]


@iet_pass
def lower_async_callables(iet, track=None, root=None, sregistry=None):
    if not isinstance(iet, CudaHostFuncCallable):
        return iet, {}

    n = len(track)

    # The `cfields` are the constant fields, that is the fields whose value
    # definitely never changes across different executions of `ìet`; the
    # `ncfields` are instead the non-constant fields, that is the fields whose
    # value may or may not change across different calls to `iet`
    indexeds = FindSymbols('indexeds|indexedbases').visit(iet)

    # Create Function -> n-dimensional array casts
    # E.g. `float (*u)[.] = (float (*)[.]) u_vec->data`
    # NOTE: a cast is needed only if the underlying data object isn't already
    # defined inside the kernel, which happens, for example, when:
    # (i) Dereferencing a PointerArray, e.g., `float (*r0)[.] = (float(*)[.]) pr0[.]`
    # (ii) Declaring a raw pointer, e.g., `float * r0 = NULL; *malloc(&(r0), ...)
    defines = set(FindSymbols('defines').visit(iet))
    bases = sorted({i.base for i in indexeds}, key=lambda i: i.name)
    casts = [PointerCast(i.function, obj=i) for i in bases
             if i not in defines]

    # need to rewrite these extra parameters into every call for this Callable
    extra_parameters = tuple([i for i in FindSymbols('basics').visit(casts)
                              if i not in iet.parameters])
    track[iet.name].extra_args = extra_parameters

    fields = iet.parameters + extra_parameters
    defines = FindSymbols('defines').visit(root.body)
    ncfields, cfields = split(fields, lambda i: i in defines)

    # SharedData -- that is the data structure that will be used by the
    # main thread to pass information down to the child thread(s)
    sdata = track[iet.name].sdata = CudaSharedData(
        name='sdata',
        npthreads=1,
        cfields=cfields,
        ncfields=ncfields,
        pname='tsdata%d' % n)
    sbase = sdata.symbolic_base

    # Prepend the SharedData fields available upon thread activation
    preactions = [
        Call("nvtxRangePush", ("__FUNCTION__",))
    ]

    # Append the flag reset
    postactions = [List(body=[
        BlankLine,
        c.Comment("Free the data block"),
        Call("free", sdata),
        Call("nvtxRangePop")
    ])]

    wrap = List(body=preactions + list(iet.body.body) + postactions)

    # pthread functions expect exactly one argument of type void*
    tparameter = Pointer(name='_%s' % sdata.name, dtype=np.void)

    # Unpack `sdata`
    unpacks = [PointerCast(sdata, tparameter), BlankLine]
    for i in flatten([cfields, ncfields]):
        if i.is_AbstractFunction:
            unpacks.append(Dereference(i, sdata))
        else:
            unpacks.append(DummyExpr(i, FieldFromPointer(i.name, sbase), init=True))

    body = iet.body._rebuild(body=[wrap, Return(Null)], unpacks=unpacks, casts=casts)
    iet = ThreadCallable(iet.name, body, tparameter)

    return iet, {}


@iet_pass
def lower_async_calls(iet, track=None, sregistry=None):
    # Definitely there won't be AsyncCalls within ThreadCallables
    if isinstance(iet, ThreadCallable):
        return iet, {}

    # Create efuncs to initialize the SharedData objects
    efuncs = OrderedDict()
    for n in FindNodes(AsyncCall).visit(iet):
        if n.name in efuncs:
            continue

        assert n.name in track
        b = track[n.name]

        sdata = b.sdata
        sbase = sdata.symbolic_base
        name = sregistry.make_name(prefix='init_%s' % sdata.name)
        body = [DummyExpr(FieldFromPointer(i._C_name, sbase), i._C_symbol)
                for i in sdata.cfields]
        parameters = sdata.cfields + (sdata,)
        efuncs[n.name] = Callable(name, body, 'void', parameters, 'static')

    # Transform AsyncCalls
    nqueues = 1  # Number of allocated asynchronous queues so far
    initialization = []
    finalization = []
    mapper = {}
    for n in FindNodes(CudaHostFuncCall).visit(iet):
        # Create `sdata` and `threads` objects for `n`
        b = track[n.name]
        name = sregistry.make_name(prefix='sdata')
        sdata = b.sdata._rebuild(name=name)
        name = sregistry.make_name(prefix='threads')

        # Call to `sdata` initialization Callable
        sbase = sdata.symbolic_base
        d = 0
        arguments = []
        for a in n.arguments + b.extra_args:
            if a in sdata.ncfields:
                continue
            elif isinstance(a, QueueID):
                # Different pthreads use different queues
                arguments.append(nqueues + d)
            else:
                arguments.append(a)
        # Each pthread has its own SharedData copy
        arguments.append(sbase + d)

        call0 = Call(efuncs[n.name].name, arguments)

        initialization.append(Definition(sdata, None, None, NullPointer()))

        # Activation
        d = 0

        activation = [c.Comment("Allocate a new block of data for this invocation"),
                      Call("posix_memalign", (VOID(Byref(sdata), '**'),
                                              64,
                                              SizeOf(sdata._C_typedata))),
                      call0]
        activation.extend([DummyExpr(FieldFromComposite(i.name, sdata[d]), i)
                           for i in sdata.ncfields])

        activation.append(
            c.Statement("cudaLaunchHostFunc(%s, (cudaHostFn_t)%s, %s)" % (
                n.stream if n.stream is not None else 0,
                n.name,
                ccode(sbase + d))),
        )

        activation = List(
            header=[c.Line(), c.Comment("Activate background task")],
            body=activation,
            footer=c.Line()
        )
        mapper[n] = activation

    if mapper:
        # Inject activation
        iet = Transformer(mapper).visit(iet)

        # Inject initialization and finalization
        initialization.append(BlankLine)
        finalization.insert(0, BlankLine)
        body = iet.body._rebuild(body=initialization + list(iet.body.body) + finalization)
        iet = iet._rebuild(body=body)
    else:
        assert not initialization
        assert not finalization

    return iet, {'efuncs': tuple(efuncs.values())}

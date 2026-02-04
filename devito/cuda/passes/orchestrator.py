from devito.symbolics import uxreplace, search
from collections import OrderedDict
from typing import Iterable

import cgen as c
import sympy

from devito.cuda.lang import CudaBB
from devito.cuda.nodes import (
    CudaChecked,
    CudaStorage,
    CudaTransferDirection,
    HostStream,
    KernelStream,
    MemCopyStream,
)
from devito.cuda.types import CudaEvent, NullPointer
from devito.ir.iet.efunc import EntryFunction
from devito.ir.iet.nodes import (
    BlankLine,
    Block,
    Call,
    Definition,
    List,
    SyncSpot,
)
from devito.ir.iet.visitors import FindNodes, Transformer, Uxreplace
from devito.ir.support.properties import AFFINE, PARALLEL
from devito.ir.support.syncs import (
    FetchUpdate,
    PrefetchUpdate,
    ReleaseLock,
    WaitLock,
    WithLock,
)
from devito.logger import debug
from devito.passes.iet.engine import iet_pass
from devito.passes.iet.orchestration import Orchestrator
from devito.symbolics.printer import ccode
from devito.tools.utils import as_mapper, filter_ordered

__all__ = ["CudaOrchestrator"]


class CudaOrchestrator(Orchestrator):
    lang = CudaBB

    _memcpy_stream = MemCopyStream()
    _host_stream = HostStream()
    _kernel_stream = KernelStream()

    def _make_waitlock(self, iet, sync_ops):
        waitloop = List(
            header=c.Comment(
                "Wait for `%s` to be copied to the host"
                % ",".join(s.function.name for s in sync_ops)
            ),
            body=[
                self.lang._map_wait_recreate_event(s, stream=self._kernel_stream)
                for s in sync_ops
            ],
            footer=c.Line(),
        )

        iet = List(body=(waitloop,) + iet.body)

        return iet, []

    def _make_withlock(self, iet, sync_ops):
        preactions = [c.Comment("Block the copy until it's safe"), BlankLine]

        # the main kernel stream should mark this as the appropriate place for it to start
        preactions.extend(
            [self.lang._map_fire_event(s, stream=self._kernel_stream) for s in sync_ops]
        )
        # these should run on the memcpy stream, so it needs to wait
        preactions.extend(
            [self.lang._map_wait_event(s, stream=self._memcpy_stream) for s in sync_ops]
        )

        # then recreate the event
        preactions.extend([self.lang._map_recreate_event(s) for s in sync_ops])

        postactions = [BlankLine, c.Comment("Raise the event")]

        # Turn `iet` into a CUDA 3D memcpy
        transfer = _make_multidimensional_async_copy(
            iet.body, CudaTransferDirection.D2H, self._memcpy_stream
        )
        efuncs = []
        postactions.extend(
            [self.lang._map_fire_event(s, stream=self._memcpy_stream) for s in sync_ops]
        )

        # The corresponding AsyncCall
        body = preactions + transfer + postactions

        iet = List(body=body)

        return iet, efuncs

    def _make_fetchupdate(self, iet, sync_ops):
        copy = _make_multidimensional_async_copy(
            iet.body, CudaTransferDirection.H2D, self._kernel_stream, dim=sync_ops[0].dim
        )

        # Perform initial fetch by the main thread
        iet = List(header=c.Comment("Initialize data stream"), body=copy)

        return iet, []

    def _make_prefetchupdate(self, iet, sync_ops):
        for op in sync_ops:
            debug(
                "copy %s to %s, size %s, dim %s, tstore %s",
                op.function,
                op.target,
                op.size,
                op.dim,
                op.tstore,
            )
        preactions = []
        preactions.extend(
            [self.lang._map_fire_event(s, stream=self._kernel_stream) for s in sync_ops]
        )
        preactions.extend(
            [self.lang._map_wait_event(s, stream=self._host_stream) for s in sync_ops]
        )
        preactions.extend([self.lang._map_recreate_event(s) for s in sync_ops])

        postactions = []
        postactions.extend(
            [self.lang._map_fire_event(s, stream=self._host_stream) for s in sync_ops]
        )
        postactions.extend(
            [
                self.lang._map_wait_recreate_event(s, stream=self._memcpy_stream)
                for s in sync_ops
            ]
        )
        postactions.extend(
            [
                self.lang._map_update_device_async(s.target, qid=self._memcpy_stream)
                for s in sync_ops
            ]
        )
        postactions.extend(
            [self.lang._map_fire_event(s, stream=self._memcpy_stream) for s in sync_ops]
        )

        # Turn `iet` into a CUDA 3D memcpy
        transfer = _make_multidimensional_async_copy(
            iet.body, CudaTransferDirection.H2D, self._memcpy_stream
        )
        efuncs = []
        preactions = []
        preactions.extend(
            [self.lang._map_fire_event(s, stream=self._kernel_stream) for s in sync_ops]
        )
        preactions.extend(
            [self.lang._map_wait_event(s, stream=self._memcpy_stream) for s in sync_ops]
        )
        preactions.extend([self.lang._map_recreate_event(s) for s in sync_ops])

        postactions = []
        postactions.extend(
            [self.lang._map_fire_event(s, stream=self._memcpy_stream) for s in sync_ops]
        )

        iet = List(body=preactions + transfer + postactions)

        return iet, efuncs

    def _nop(self, _iet, _sync_ops):
        return _iet, []

    def _replace_locks(self, iet):
        # replace locks with CUDA events
        lock_mapper = {}

        for n in FindNodes(SyncSpot).visit(iet):
            replace_ops = []
            for s in n.sync_ops:
                if s.handle:
                    s.handle = CudaEvent(s.handle.name)
                replace_ops.append(s)

            lock_mapper[n] = SyncSpot(replace_ops, n.body)

        iet = Uxreplace(lock_mapper).visit(iet)

        return iet

    @iet_pass
    def process(self, iet):
        iet = self._replace_locks(iet)

        sync_spots = FindNodes(SyncSpot).visit(iet)

        check_exception_raised = [c.Line("devicerm |= exceptionOccured();")]

        if not sync_spots:
            if isinstance(iet, EntryFunction):
                # yuck
                iet = iet._rebuild(
                    body=List(
                        body=[
                            iet.body,
                            CudaChecked(Call("cudaStreamSynchronize", KernelStream())),
                            CudaChecked(Call("cudaStreamSynchronize", MemCopyStream())),
                            CudaChecked(Call("cudaStreamSynchronize", HostStream())),
                        ]
                        + check_exception_raised
                    )
                )
            return iet, {}

        callbacks = OrderedDict(
            [
                (WithLock, self._make_withlock),
                (WaitLock, self._make_waitlock),
                (FetchUpdate, self._make_fetchupdate),
                (PrefetchUpdate, self._make_prefetchupdate),
                (ReleaseLock, self._nop),
            ]
        )

        # The SyncOps are to be processed in a given order
        key = lambda s: list(callbacks).index(s)

        efuncs = []
        subs = {}
        events = []
        for n in sync_spots:
            mapper = as_mapper(n.sync_ops, lambda i: type(i))
            for t in sorted(mapper, key=key):
                events.extend([e.handle for e in mapper[t] if e.handle is not None])
                subs[n], v = callbacks[t](subs.get(n, n), mapper[t])
                efuncs.extend(v)

        iet = Transformer(subs).visit(iet)

        events = [
            List(
                body=[
                    Definition(e, None, None, NullPointer()),
                    self.lang.mapper["create-event"](e._C_symbol),
                ]
            )
            for e in filter_ordered(events)
        ]
        iet = iet._rebuild(
            body=List(
                body=[
                    events,
                    iet.body,
                    CudaChecked(Call("cudaStreamSynchronize", KernelStream())),
                    CudaChecked(Call("cudaStreamSynchronize", MemCopyStream())),
                    CudaChecked(Call("cudaStreamSynchronize", HostStream())),
                ]
                + check_exception_raised
            )
        )

        return iet, {"efuncs": efuncs}


def is_multidimensional_memcpy(iet) -> bool:
    if isinstance(iet, Iterable):
        return all(is_multidimensional_memcpy(n) for n in iet)

    if iet.is_Iteration:
        if (
            iet.dim.is_Modulo
            or not all(p in iet.properties for p in [AFFINE, PARALLEL])
            or not is_multidimensional_memcpy(iet.children)
        ):
            return False
    elif iet.is_ExpressionBundle:
        if not all(is_memcpy_relaxed(e) for e in iet.exprs):
            return False
    elif iet.is_Expression:
        if not is_memcpy_relaxed(iet):
            return False
    else:
        return False

    return True


def is_memcpy_relaxed(e) -> bool:
    return (
        (e.expr.lhs.function.is_Array or e.expr.lhs.function.is_DiscreteFunction)
        and e.expr.lhs.function._mem_mapped
        and e.expr.rhs.is_Indexed
        and not e.expr.lhs.function.is_SparseFunction
        and not e.expr.rhs.function.is_SparseFunction
    )


def gather_multidimensional_memcpy(iet) -> list[tuple]:
    dims = list()
    exprs = list()

    if isinstance(iet, Iterable):
        for n in iet:
            sub_dims, sub_exprs = gather_multidimensional_memcpy(n)
            dims.extend(sub_dims)
            exprs.extend(sub_exprs)

        return dims, exprs

    if iet.is_Iteration and all(p in iet.properties for p in [PARALLEL]):
        dims.append((iet.dimensions, iet.limits))
        sub_dims, sub_exprs = gather_multidimensional_memcpy(iet.children)
        dims.extend(sub_dims)
        exprs.extend(sub_exprs)
    elif iet.is_ExpressionBundle:
        exprs.extend([e.expr for e in iet.exprs if is_memcpy_relaxed(e)])
    elif iet.is_Expression and is_memcpy_relaxed(iet):
        exprs.extend(iet.expr)

    return dims, exprs


# Search for expressions of the form "x % 1", which can be replaced with a literal 0
def q_mod_can_elide(expr):
    if isinstance(expr, sympy.Mod) and expr.args[1] == 1:
        return True


def _make_multidimensional_async_copy(iet, direction, stream, dim=None) -> list[Block]:
    dims, exprs = gather_multidimensional_memcpy(iet)

    ret = []
    debug("examined iet %s" % str(iet))
    debug("found %d expressions for multidimensional async copy" % len(exprs))

    mapper = {}
    if dim is not None:
        # if this dimension has zero size, replace it with zero in all expressions
        if dim.symbolic_min == dim.symbolic_max:
            mapper = {dim: sympy.Integer(0)}

    mods = search(exprs, q_mod_can_elide, "all", "dfs")
    for mod in mods:
        mapper[mod] = sympy.Integer(0)
        debug("eliding modulo operation %s" % str(mod))

    for k, v in mapper.items():
        debug("replacing %s with %s" % (str(k), str(v)))
    if len(mapper) > 0:
        new_exprs = []
        for expr in exprs:
            new_exprs.append(uxreplace(expr, mapper))

        exprs = new_exprs

    debug("found %d expressions for multidimensional async copy" % len(exprs))
    for expr in exprs:
        src = expr.rhs
        dst = expr.lhs

        src_storage = CudaStorage(src.function)
        dst_storage = CudaStorage(dst.function)

        call = Call(
            "devito::cuda::async_buffer_copy_slice",
            [
                dst_storage.function._C_name,
                src_storage.function._C_name,
                dst.indices[-4] if len(dst.indices) >= 4 else 0,
                src.indices[-4] if len(src.indices) >= 4 else 0,
                "cudaMemcpyHostToDevice"
                if direction == CudaTransferDirection.H2D
                else "cudaMemcpyDeviceToHost",
                ccode(stream),
            ],
        )

        ret.append(call)

    return ret

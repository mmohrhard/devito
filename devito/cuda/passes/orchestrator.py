from typing import Iterable
import cgen as c
import sympy

from collections import OrderedDict

from devito.cuda.types import CudaEvent, NullPointer
from devito.ir.iet.efunc import EntryFunction
from devito.ir.iet.nodes import (
    BlankLine,
    Block,
    Call,
    Callable,
    Definition,
    List,
    SyncSpot,
)
from devito.ir.iet.utils import derive_parameters
from devito.ir.iet.visitors import FindNodes, Transformer, Uxreplace
from devito.ir.support.properties import AFFINE, PARALLEL
from devito.ir.support.syncs import (
    FetchUpdate,
    PrefetchUpdate,
    ReleaseLock,
    WaitLock,
    WithLock,
)
from devito.passes.iet.engine import iet_pass
from devito.passes.iet.orchestration import Orchestrator
from devito.symbolics.printer import ccode
from devito.tools.utils import as_mapper, as_tuple, filter_ordered, flatten

from devito.cuda.lang import CudaBB
from devito.cuda.nodes import (
    CudaChecked,
    CudaStorage,
    CudaTransferDirection,
    HostStream,
    KernelStream,
    MemCopyStream,
)

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
        postactions = [self.lang._map_update_device(s.target, s.imask) for s in sync_ops]

        # Turn init IET into a Callable
        name = self.sregistry.make_name(prefix="init_device")
        body = List(body=iet.body + tuple(postactions))
        parameters = derive_parameters(body)
        efunc = Callable(name, body, "void", parameters, "static")

        # Perform initial fetch by the main thread
        iet = List(
            header=c.Comment("Initialize data stream"), body=Call(name, parameters)
        )

        return iet, [efunc]

    def _make_prefetchupdate(self, iet, sync_ops):
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

    if iet.is_Iteration and all(p in iet.properties for p in [AFFINE, PARALLEL]):
        dims.append((iet.dimensions, iet.limits))
        sub_dims, sub_exprs = gather_multidimensional_memcpy(iet.children)
        dims.extend(sub_dims)
        exprs.extend(sub_exprs)
    elif iet.is_ExpressionBundle:
        exprs.extend([e.expr for e in iet.exprs if is_memcpy_relaxed(e)])
    elif iet.is_Expression and is_memcpy_relaxed(iet):
        exprs.extend(iet.expr)

    return dims, exprs


def _make_multidimensional_async_copy(iet, direction, stream) -> list[Block]:
    dims, exprs = gather_multidimensional_memcpy(iet)

    ret = []

    for expr in exprs:
        src = expr.rhs
        dst = expr.lhs

        d = flatten(x[0] for x in dims)
        iteration_dimensions = list(x for x in d if not x.is_Derived)
        iterators = list(x for x in d if x.is_Derived)

        src_non_iterated_dimensions = [
            x
            for x in src.function.dimensions
            if x not in iteration_dimensions and x not in iterators
        ]
        dst_non_iterated_dimensions = [
            x
            for x in dst.function.dimensions
            if x not in iteration_dimensions and x not in iterators
        ]

        src_storage = CudaStorage(src.function)
        dst_storage = CudaStorage(dst.function)

        src_strides = [
            sympy.Mul(
                *(
                    src.function.symbolic_shape[i + 1 :]
                    if i < len(src.indices) - 1
                    else [1]
                )
            )
            for i in range(len(src.indices))
        ]
        src_offset = sympy.Add(
            *[
                sympy.Mul(src_strides[i], src.indices[i])
                for i in range(len(src.indices))
                if src.function.dimensions[i] in src_non_iterated_dimensions
            ]
        )
        src_ptr = (
            src_storage.host_storage
            if direction == CudaTransferDirection.H2D
            else src_storage.device_storage
        )

        src_ptr = "&((float *)(%s))[%s]" % (src_ptr, ccode(src_offset))
        src_pitch = "sizeof(float) * (%s)" % ccode(src.function.symbolic_shape[-1])
        src_pos = "make_cudaPos(sizeof(float) * (%s), (%s), (%s))" % (
            as_tuple(
                [
                    (
                        ccode(
                            dims[-i][1][0]
                            + (
                                (src.indices[-i] - iterators[-i])
                                if len(iterators) > (i - 1)
                                else 0
                            )
                            - src.function.dimensions[-i].symbolic_min
                        )
                        if len(dims) >= i
                        else "0"
                    )
                    for i in range(1, 4)
                ]
            )
        )

        src_x = ccode(src.function.symbolic_shape[-1])
        src_y = ccode(src.function.symbolic_shape[-2])

        dst_ptr = (
            dst_storage.device_storage
            if direction == CudaTransferDirection.H2D
            else dst_storage.host_storage
        )

        dst_strides = [
            sympy.Mul(
                *(
                    dst.function.symbolic_shape[i + 1 :]
                    if i < len(dst.indices) - 1
                    else [1]
                )
            )
            for i in range(len(dst.indices))
        ]
        dst_offset = sympy.Add(
            *[
                sympy.Mul(dst_strides[i], dst.indices[i])
                for i in range(len(dst.indices))
                if dst.function.dimensions[i] in dst_non_iterated_dimensions
            ]
        )

        dst_ptr = "&((float *)(%s))[%s]" % (
            dst_ptr,
            ccode(dst_offset),
        )
        dst_pitch = "sizeof(float) * (%s)" % ccode(dst.function.symbolic_shape[-1])
        dst_pos = "make_cudaPos(sizeof(float) * (%s), (%s), (%s))" % (
            as_tuple(
                [
                    (
                        ccode(
                            dims[-i][1][0]
                            + (
                                (dst.indices[-i] - iterators[-i])
                                if len(iterators) > (i - 1)
                                else 0
                            )
                            - dst.function.dimensions[-i].symbolic_min
                        )
                        if len(dims) >= i
                        else "0"
                    )
                    for i in range(1, 4)
                ]
            )
        )

        dst_x = ccode(dst.function.symbolic_shape[-1])
        dst_y = ccode(dst.function.symbolic_shape[-2])

        assert len(dims) <= max(
            len(src.function.dimensions), len(dst.function.dimensions)
        )

        extent = "make_cudaExtent(sizeof(float) * (%s), %s, %s)" % (
            ccode(dims[-1][1][1] - dims[-1][1][0]),
            (
                ("(" + ccode(dims[-2][1][1] - dims[-2][1][0]) + ")")
                if len(dims) >= 2
                else "1"
            ),
            (
                ("(" + ccode(dims[-3][1][1] - dims[-3][1][0]) + ")")
                if len(dims) >= 3
                else "1"
            ),
        )

        ops = [
            c.Statement("struct cudaMemcpy3DParms copy_params = {0}"),
            c.Statement(
                "copy_params.srcPtr = make_cudaPitchedPtr(%s, %s, %s, %s)"
                % (src_ptr, src_pitch, src_x, src_y)
            ),
            c.Statement("copy_params.srcPos = %s" % src_pos),
            c.Statement(
                "copy_params.dstPtr = make_cudaPitchedPtr(%s, %s, %s, %s)"
                % (dst_ptr, dst_pitch, dst_x, dst_y)
            ),
            c.Statement("copy_params.dstPos = %s" % dst_pos),
            c.Statement("copy_params.extent = %s" % extent),
            c.Statement("copy_params.kind = cudaMemcpyDefault"),
            c.Statement(
                "CudaChecked(cudaMemcpy3DAsync(&copy_params, %s))" % ccode(stream)
            ),
        ]

        ret.append(Block(body=ops))

    return ret

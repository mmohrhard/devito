import cgen as c

from collections import OrderedDict

from devito.cuda.types import CudaEvent, NullPointer
from devito.cuda.utils import cuda_derive_parameters
from devito.ir.iet.efunc import EntryFunction
from devito.ir.iet.nodes import BlankLine, Call, Callable, Definition, List, SyncSpot
from devito.ir.iet.utils import derive_parameters
from devito.ir.iet.visitors import FindNodes, Transformer, Uxreplace
from devito.ir.support.syncs import (
    FetchUpdate,
    PrefetchUpdate,
    ReleaseLock,
    WaitLock,
    WithLock,
)
from devito.passes.iet.engine import iet_pass
from devito.passes.iet.orchestration import Orchestrator
from devito.tools.utils import as_mapper, filter_ordered

from devito.cuda.lang import CudaBB
from devito.cuda.nodes import (
    CudaChecked,
    CudaHostFuncCall,
    CudaHostFuncCallable,
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
        preactions.extend(
            [
                self.lang._map_update_host_async(s.function, qid=self._memcpy_stream)
                for s in sync_ops
            ]
        )
        preactions.extend(
            [self.lang._map_fire_event(s, stream=self._memcpy_stream) for s in sync_ops]
        )

        # these should run on the host stream (not default stream)
        preactions.extend(
            [
                self.lang._map_wait_recreate_event(s, stream=self._host_stream)
                for s in sync_ops
            ]
        )
        postactions = [BlankLine, c.Comment("Raise the event")]
        postactions.extend(
            [self.lang._map_fire_event(s, stream=self._host_stream) for s in sync_ops]
        )

        # Turn `iet` into an AsyncCallable so that subsequent passes know
        # that we're happy for this Callable to be executed asynchronously
        name = self.sregistry.make_name(prefix="copy_device_to_host")
        async_body = List(body=iet.body)
        parameters = cuda_derive_parameters(async_body)
        async_body = async_body._rebuild()
        efunc = CudaHostFuncCallable(name, async_body, parameters=parameters)

        # The corresponding AsyncCall
        body = (
            preactions
            + [CudaHostFuncCall(name, efunc.parameters, stream=self._host_stream)]
            + postactions
        )

        iet = List(body=body)

        return iet, [efunc]

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

        # Turn `iet` into an AsyncCallable so that subsequent passes know
        # that we're happy for this Callable to be executed asynchronously
        name = self.sregistry.make_name(prefix="prefetch_host_to_device")
        body = iet.body

        parameters = cuda_derive_parameters(body)
        efunc = CudaHostFuncCallable(name, body, parameters=parameters)

        # The corresponding AsyncCall
        iet = List(
            body=preactions
            + [CudaHostFuncCall(name, efunc.parameters, stream=self._host_stream)]
            + postactions
        )

        return iet, [efunc]

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

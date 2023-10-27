from functools import singledispatch
from devito.types.dense import AliasFunction
from devito.types.parallel import DeviceCreate, DeviceRM, UpdateDevice, UpdateHost

from devito.tools.utils import filter_sorted, flatten

from devito.symbolics.printer import ccode
from devito.symbolics.extended_sympy import (
    VOID,
    Byref,
    CondOr,
    IndexedPointer,
    Keyword,
    ReservedWord,
    SizeOf,
)

from devito.ir.iet.nodes import Call, Conditional, Definition, List, ListInitializer
from devito.ir.iet.efunc import EntryFunction
from devito.ir.iet.visitors import FindNodes, FindSymbols, MapExprStmts, Transformer

from devito.passes.iet import DataManager
from devito.passes.iet.definitions import Storage
from devito.passes.iet.parpragma import PragmaTransfer
from devito.passes.iet.engine import iet_pass, iet_visit
from devito.passes.iet.misc import is_on_device

from devito.cuda.passes.tuning import kernel_tuning
from devito.cuda.utils import flatten_dict
from devito.cuda.lang import CudaBB
from devito.cuda.nodes import (
    CudaCall,
    CudaCallable,
    CudaDealloc,
    CudaKernelPointerCast,
    CudaTransfer,
    DeviceCall,
    DeviceFunction,
)

__all__ = ["DeviceCudaDataManager"]


class DeviceCudaDataManager(DataManager):
    lang = CudaBB

    def __init__(self, sregistry, options):
        """
        Parameters
        ----------
        sregistry : SymbolRegistry
            The symbol registry, to quickly access the special symbols that may
            appear in the IET.
        options : dict
            The optimization options.
            Accepted: ['gpu-fit'].
            * 'gpu-fit': an iterable of `Function`s that are guaranteed to fit
              in the device memory. By default, all `Function`s except saved
              `TimeFunction`'s are assumed to fit in the device memory.
        """
        super().__init__(sregistry)
        self.gpu_fit = options["gpu-fit"]
        self.gpu_nofit = options["gpu-nofit"]

    def _alloc_local_array_on_high_bw_mem(self, site, obj, storage, devicerm=None):
        """
        Allocate a local Array in the device high bandwidth memory.
        """
        nbytes = SizeOf(obj._C_typedata) * obj.size

        allocs = Call(
            "PER_DEVICE_TEMP_GET",
            (ReservedWord(str(obj._C_typedata)), obj._C_symbol, nbytes),
        )

        free = Call("PER_DEVICE_TEMP_DESTROY", (obj._C_name,))
        free = Conditional(DeviceRM(), free)

        storage.update(obj, site, allocs=allocs, frees=free)

    def _alloc_object_on_low_lat_mem(self, site, obj, storage):
        """
        Allocate a LocalObject in the low latency memory.
        """
        decl = Definition(obj, cargs=obj.cargs)

        if obj._C_init:
            definition = (decl, obj._C_init)
        else:
            definition = decl

        frees = obj._C_free

        storage.update(obj, site, objs=definition, frees=frees)

    def _alloc_array_on_low_lat_mem(self, site, obj, storage):
        """
        Allocate an Array in the low latency memory.
        """
        shape = "".join("[%s]" % ccode(i) for i in obj.symbolic_shape)
        alignment = self.lang["aligned"](obj._data_alignment)
        if obj.initvalue is None:
            initvalue = None
        else:
            initvalue = ListInitializer(obj.initvalue)
        alloc = Definition(obj, shape=shape, qualifier=alignment, initvalue=initvalue)

        storage.update(obj, site, allocs=alloc)

    def _alloc_scalar_on_low_lat_mem(self, site, expr, storage):
        """
        Allocate a Scalar in the low latency memory.
        """
        storage.map(expr.write, site, expr, expr._rebuild(init=True))

    def _alloc_host_array_on_high_bw_mem(self, site, obj, storage, *args):
        """
        Allocate a host Array in the host high bandwidth memory.
        """
        decl = Definition(obj)

        memptr = VOID(Byref(obj._C_symbol), "**")
        alignment = obj._data_alignment
        nbytes = SizeOf(obj._C_typedata) * obj.size
        alloc = self.lang["host-alloc"](memptr, alignment, nbytes)

        free = self.lang["host-free"](obj._C_symbol)

        storage.update(obj, site, allocs=(decl, alloc), frees=free)

    def _alloc_mapped_array_on_high_bw_mem(self, site, obj, storage, *args):
        """
        Allocate a mapped Array in the host high bandwidth memory.
        """

        nbytes_arg = SizeOf(obj.indexed._C_typedata) * obj.size

        alloc = List(
            body=[
                Call(
                    "PER_DEVICE_ARRAY_TEMP_DECLARE",
                    (obj._C_symbol, ReservedWord(obj._C_typedata)),
                ),
                Call(
                    "PER_DEVICE_ARRAY_TEMP_GET", (obj._C_symbol, nbytes_arg), retobj=obj
                ),
            ]
        )

        free = Conditional(
            DeviceRM(), Call("PER_DEVICE_ARRAY_TEMP_DESTROY", (obj._C_symbol,))
        )

        storage.update(obj, site, allocs=alloc, frees=free)

    def _alloc_object_array_on_low_lat_mem(self, site, obj, storage):
        """
        Allocate an Array of Objects in the low latency memory.
        """
        shape = "".join("[%s]" % ccode(i) for i in obj.symbolic_shape)
        decl = Definition(obj, shape=shape)

        storage.update(obj, site, allocs=decl)

    def _alloc_pointed_array_on_high_bw_mem(self, site, obj, storage):
        """
        Allocate the following objects in the high bandwidth memory:

            * The pointer array `obj`;
            * The pointee Array `obj.array`

        If the pointer array is defined over `sregistry.threadid`, that is a thread
        Dimension, then each `obj.array` slice is allocated and freed individually
        by the owner thread.
        """
        # The pointer array
        decl = Definition(obj)

        memptr = VOID(Byref(obj._C_symbol), "**")
        alignment = obj._data_alignment
        nbytes = SizeOf(Keyword("%s*" % obj._C_typedata)) * obj.dim.symbolic_size
        alloc0 = self.lang["host-alloc"](memptr, alignment, nbytes)

        free0 = self.lang["host-free"](obj._C_symbol)

        # The pointee Array
        pobj = IndexedPointer(obj._C_symbol, obj.dim)
        memptr = VOID(Byref(pobj), "**")
        nbytes = SizeOf(obj._C_typedata) * obj.array.size
        alloc1 = self.lang["host-alloc"](memptr, alignment, nbytes)

        free1 = self.lang["host-free"](pobj)

        # Dump
        if obj.dim is self.sregistry.threadid:
            storage.update(
                obj,
                site,
                allocs=(decl, alloc0),
                frees=free0,
                pallocs=(obj.dim, alloc1),
                pfrees=(obj.dim, free1),
            )
        else:
            storage.update(obj, site, allocs=(decl, alloc0, alloc1), frees=(free0, free1))

    def _map_array_on_high_bw_mem(self, _site, _obj, _storage):
        """
        Map an Array already defined in the host memory in to the device high
        bandwidth memory.
        """

        # When using CUDA we allocate everything in a device-visible manner
        return

    def _map_function_on_high_bw_mem(
        self,
        site,
        obj,
        storage,
        devicerm,
        read_only=False,
        devicecreate=None,
        updatehost=None,
        updatedevice=None,
    ):
        """
        Map a Function already defined in the host memory in to the device high
        bandwidth memory.

        Notes
        -----
        In essence, the difference between `_map_function_on_high_bw_mem` and
        `_map_array_on_high_bw_mem` is that the former triggers a data transfer to
        synchronize the host and device copies, while the latter does not.
        """
        if devicecreate:
            mmap = [
                self.lang._map_update_device(
                    obj, condition=CondOr(devicecreate, updatedevice)
                )
            ]
        else:
            mmap = self.lang._map_to(obj)

        if read_only is False:
            unmap = [
                self.lang._map_update_host(obj, condition=CondOr(updatehost, devicerm)),
                self.lang._map_release(obj, devicerm=devicerm),
            ]
        else:
            unmap = self.lang._map_delete(obj, devicerm=devicerm)

        storage.update(obj, site, maps=mmap, unmaps=unmap)

    def _dump_transfers(self, iet, storage):
        mapper = {}
        for k, v in storage.items():
            if v.maps or v.unmaps:
                mapper[iet.body] = iet.body._rebuild(
                    maps=flatten(v.maps), unmaps=flatten(v.unmaps)
                )

        processed = Transformer(mapper, nested=True).visit(iet)

        return processed

    @iet_visit
    def derive_transfers(self, iet):
        """
        Collect all symbols that cause host-device data transfer, distinguishing
        between reads and writes.
        """

        def needs_transfer(f):
            return (
                f._mem_mapped
                and not isinstance(f, AliasFunction)
                and is_on_device(f, self.gpu_fit)
            )

        writes = set()
        reads = set()
        for i, v in MapExprStmts().visit(iet).items():
            if (
                not any(isinstance(j, self.lang.DeviceIteration) for j in v)
                and not isinstance(i, DeviceCall)
                and not isinstance(iet, DeviceFunction)
            ):
                # Not an offloaded Iteration tree
                continue

            writes.update({w for w in i.writes if needs_transfer(w)})
            reads.update(
                {f for f in i.functions if needs_transfer(f) and f not in writes}
            )

        return (reads, writes)

    @iet_pass
    def place_transfers(self, iet, **kwargs):
        """
        Create a new IET with host-device data transfers. This requires mapping
        symbols to the suitable memory spaces.
        """

        @singledispatch
        def _place_transfers(iet, mapper):
            return iet, {}

        @_place_transfers.register(EntryFunction)
        def _(iet, mapper):
            try:
                reads, writes = list(zip(*mapper.values()))
            except ValueError:
                return iet, {}
            reads = set(flatten(reads))
            writes = set(flatten(writes))

            # Special symbol which gives user code control over data deallocations
            devicecreate = DeviceCreate()
            devicerm = DeviceRM()
            updatehost = UpdateHost()
            updatedevice = UpdateDevice()
            storage = Storage()
            for i in filter_sorted(writes):
                if i.is_Array:
                    self._map_array_on_high_bw_mem(iet, i, storage)
                else:
                    self._map_function_on_high_bw_mem(
                        iet,
                        i,
                        storage,
                        devicerm,
                        devicecreate=devicecreate,
                        updatehost=updatehost,
                        updatedevice=updatedevice,
                    )
            for i in filter_sorted(reads - writes):
                if i.is_Array:
                    self._map_array_on_high_bw_mem(iet, i, storage)
                else:
                    self._map_function_on_high_bw_mem(
                        iet,
                        i,
                        storage,
                        devicerm,
                        read_only=True,
                        devicecreate=devicecreate,
                        updatehost=updatehost,
                        updatedevice=updatedevice,
                    )

            iet = self._dump_transfers(iet, storage)

            return iet, {}

        return _place_transfers(iet, mapper=kwargs["mapper"])

    @iet_visit
    def derive_cuda_casts(self, iet, **kwargs):
        # Don't generate unnecessary casts in CUDA kernels
        kernels = FindNodes(CudaCallable).visit(iet)
        mapper = {}

        for kernel in kernels:
            indexeds = FindSymbols("indexeds|indexedbases").visit(kernel)
            defines = set(FindSymbols("defines").visit(kernel)) - set(kernel.parameters)
            bases = sorted({i.base for i in indexeds}, key=lambda i: i.name)
            casts = [
                CudaKernelPointerCast(i.function, obj=i)
                for i in bases
                if i.function not in defines
            ]

            # Incorporate the newly created casts
            if casts:
                mapper[kernel] = kernel._rebuild(body=kernel.body._rebuild(casts=casts))

        return mapper

    @iet_visit
    def derive_cuda_kernel_call_parameters(self, iet, mapper: dict):
        calls = FindNodes(CudaCall).visit(iet)
        cmapper = {}
        for kernel, replacement in mapper.items():
            template_args = replacement.template_arguments
            our_calls = [c for c in calls if c.name == kernel.name]
            for c in our_calls:
                cmapper[c] = c._rebuild(template_arguments=template_args)

        return {**mapper, **cmapper}

    @iet_pass
    def place_cuda_non_kernel_casts(self, iet, **kwargs):
        if not isinstance(iet, CudaCallable):
            cuda_filter = lambda n: isinstance(
                n,
                (
                    CudaCall,
                    CudaCallable,
                    CudaDealloc,
                    PragmaTransfer,
                    CudaTransfer,
                ),
            )
            # Candidates
            indexeds = FindSymbols(
                "indexeds|indexedbases", stop_filter=cuda_filter
            ).visit(iet)

            # Create Function -> n-dimensional array casts
            # E.g. `float (*u)[.] = (float (*)[.]) u_vec->data`
            # NOTE: a cast is needed only if the underlying data object isn't already
            # defined inside the kernel, which happens, for example, when:
            # (i) Dereferencing a PointerArray,
            #     e.g., `float (*r0)[.] = (float(*)[.]) pr0[.]`
            # (ii) Declaring a raw pointer, e.g., `float * r0 = NULL; *malloc(&(r0), ...)
            # we use iet.body here because we manually futz with the defines for some
            # nodes to coerce Devito into outputting function signatures the way
            # we want them
            defines = set(
                [
                    x.function
                    for x in FindSymbols("defines", stop_filter=cuda_filter).visit(
                        iet.body
                    )
                ]
            )
            bases = sorted({i.base for i in indexeds}, key=lambda i: i.name)
            casts = [
                self.lang.PointerCast(i.function, obj=i)
                for i in bases
                if i.function not in defines
            ]

            # Incorporate the newly created casts
            if casts:
                iet = iet._rebuild(body=iet.body._rebuild(casts=casts))

        return iet, {}

    @iet_pass
    def place_cuda_casts(self, iet, **kwargs):
        return Transformer(kwargs["mapper"]).visit(iet), {}

    def process(self, graph):
        """
        Apply the `place_transfers`, `place_definitions` and `place_casts` passes.
        """
        mapper = self.derive_transfers(graph)
        self.place_transfers(graph, mapper=mapper)
        self.place_definitions(graph)
        cast_mapper = self.derive_cuda_casts(graph)
        cast_mapper = flatten_dict(cast_mapper, prefix=None)
        cast_mapper = flatten_dict(
            self.derive_cuda_kernel_call_parameters(graph, mapper=cast_mapper),
            prefix=None,
        )
        self.place_cuda_casts(graph, mapper=cast_mapper)
        self.place_cuda_non_kernel_casts(graph)
        kernel_tuning(graph)

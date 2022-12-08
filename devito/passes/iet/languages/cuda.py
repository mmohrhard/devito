from collections import OrderedDict
from ctypes import POINTER
from enum import Enum
from functools import cached_property, singledispatch
import cgen as c
from sympy import Or
from devito.ir.iet.visitors import Visitor

from devito.symbolics.extended_sympy import FieldFromComposite, Null, ReservedWord
from devito.types.misc import Pointer
from devito.tools.data_structures import Bunch, DefaultOrderedDict
from devito.ir.support.syncs import CudaFetchUpdate, CudaFireEvent, CudaPrefetchUpdate, CudaWaitEvent, CudaWithEvent, FetchUpdate, PrefetchUpdate, ReleaseLock, WaitLock, WithLock
from devito.ir.iet.efunc import AsyncCall, AsyncCallable, ThreadCallable
from devito.ir.iet.cuda import CudaTransferDirection
from devito.passes.iet.definitions import DeviceAwareDataManager, Storage
from devito.passes.iet.langbase import make_sections_from_imask
from devito.symbolics.printer import ccode
from devito.ir.iet.nodes import BlankLine, BusyWait, Conditional, Dereference, PointerCast, Return, Section, SyncSpot, Transfer, While
from devito.tools.utils import as_list, as_mapper, as_tuple, dtype_to_ctype, filter_sorted, flatten, split
from devito.types.parallel import CudaStream, DeviceRM, QueueID, SharedData, ThreadArray
import numpy as np
from sympy import Max
from devito.logger import warning, debug
from devito.arch import CUDA, NVIDIAX
from devito.ir import (Call, Callable, CudaCall, DeviceCall, DummyExpr, DPtr, EntryFunction, List, CudaCallable,
                       Block, ParallelIteration, ParallelTree, Pragma, Definition, Iteration, Node,
                       FindNodes, FindSymbols, Uxreplace, Transformer, Lambda, AddressOf, CLiteral,
                       MapExprStmts, DeviceFunction, make_callable, Expression, OpInc)
from devito.passes.iet.engine import iet_pass, iet_visit
from devito.passes.iet.definitions import DataManager
from devito.passes.iet.orchestration import Orchestrator
from devito.passes.iet.parpragma import (PragmaDeviceAwareTransformer, PragmaLangBB,
                                         PragmaTransfer, PragmaDeviceAwareDataManager)
from devito.symbolics import (Byref, DefFunction, FieldFromPointer, IndexedPointer,
                              ListInitializer, SizeOf, VOID, Keyword, ccode,
                               CondEq, CondNe, CondOr)
from devito.passes.iet.languages.C import CBB
from devito.passes.iet.languages.openmp import OmpRegion, OmpIteration
from devito.passes.iet.languages.utils import make_clause_reduction
from devito.passes.iet.misc import is_on_device
from devito.ir.iet.utils import derive_parameters, retrieve_iteration_tree, filter_iterations
from devito.symbolics import Macro, cast_mapper

from devito.tools import filter_ordered
from devito.types import DevicePointer, Symbol, Constant, DeviceRM, DeviceCreate, UpdateDevice, UpdateHost
from devito.types.dense import AliasFunction

__all__ = ['DeviceCudaizer', 'DeviceCudaDataManager', 'CudaOrchestrator', 'cuda_eventify']


class DeviceCudaIteration(ParallelIteration):

    @classmethod
    def _make_construct(cls, **kwargs):
        return 'acc parallel loop'

    @classmethod
    def _make_clauses(cls, ncollapse=None, reduction=None, tile=None, **kwargs):
        clauses = []

        if ncollapse:
            clauses.append('collapse(%d)' % (ncollapse or 1))
        elif tile:
            clauses.append('tile(%s)' % ','.join(str(i) for i in tile))

        if reduction:
            clauses.append(make_clause_reduction(reduction))

        indexeds = FindSymbols('indexeds').visit(kwargs['nodes'])
        deviceptrs = filter_ordered(i.name for i in indexeds if i.function._mem_local)
        presents = filter_ordered(i.name for i in indexeds
                                  if (is_on_device(i, kwargs['gpu_fit']) and
                                      i.name not in deviceptrs))

        # The NVC 20.7 and 20.9 compilers have a bug which triggers data movement for
        # indirectly indexed arrays (e.g., a[b[i]]) unless a present clause is used
        if presents:
            clauses.append("present(%s)" % ",".join(presents))

        if deviceptrs:
            clauses.append("deviceptr(%s)" % ",".join(deviceptrs))

        return clauses

    @classmethod
    def _process_kwargs(cls, **kwargs):
        kwargs = super()._process_kwargs(**kwargs)

        kwargs.pop('gpu_fit', None)

        kwargs.pop('schedule', None)
        kwargs.pop('parallel', None)
        kwargs.pop('chunk_size', None)
        kwargs.pop('nthreads', None)
        kwargs.pop('tile', None)

        return kwargs


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
        return "%s->%s" % (self.function._C_name, self.function._C_field_operator_allocated)
    
    @cached_property
    def size(self):
        return ('sizeof(%s) * ' % (self.function.indexed._C_typedata)) + '*'.join("(" + ccode(j) + ")" for i, j in self.sections)

    @property
    def imask(self):
        return self._imask

    @cached_property
    def sections(self):
        return make_sections_from_imask(self.function, self.imask)


class CudaTransfer(CudaStorage, Transfer, Node):
    """
    A data transfer between host and CUDA device.
    """

    def __init__(self, function, imask=None, condition=None, direction=CudaTransferDirection.H2D, delete=None, stream=None):
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
        super().__init__('if (cudaPeekAtLastError() != 0 ) { cudaError_t err = cudaGetLastError(); printf("\\n!E %s: %s\\n",cudaGetErrorName(err), cudaGetErrorString(err));}')

class CudaChecked(Call):
    def __init__(self, arguments=None):
        super().__init__("CudaChecked", arguments=[arguments])

    @cached_property
    def expr_symbols(self):
        return flatten([x.expr_symbols for x in flatten(self.arguments)])

class NullPointer(ReservedWord):
    def __new__(cls):
        return super().__new__(cls, "nullptr")

class CudaBB(PragmaLangBB):

    mapper = {
        # Misc
        'name': 'CUDA',
        'headers': ['cuda.h', 'cuda_runtime_api.h', 'nvtx3/nvToolsExt.h', 'stdio.h', 'assert.h', 'devito/devito_cuda.cuh'],
        # Platform mapping
        CUDA: None,
        NVIDIAX: None,
        # Runtime library
        'aligned': lambda i:
            '__attribute__((aligned(%d)))' % i,
        'init': lambda args:
            List(body=[Definition(HostStream(), initvalue=NullPointer(), prefix="static"),
                       Definition(MemCopyStream(), initvalue=NullPointer(), prefix="static"),
                       Definition(KernelStream(), initvalue=NullPointer(), prefix="static"),

                       Conditional(CondEq(HostStream(), NullPointer()), Call("cudaStreamCreateWithFlags", (Byref(HostStream()), "cudaStreamNonBlocking"))),
                       Conditional(CondEq(MemCopyStream(), NullPointer()), Call("cudaStreamCreateWithFlags", (Byref(MemCopyStream()), "cudaStreamNonBlocking"))),
                       Conditional(CondEq(KernelStream(), NullPointer()), Call("cudaStreamCreateWithFlags", (Byref(KernelStream()), "cudaStreamNonBlocking"))),

                       Call("nvtxRangePush", ("__FUNCTION__", )),
                       c.Statement(f'printf("devicerm=%d, updatehost=%d, updatedevice=%d, devicecreate=%d\\n", devicerm, updatehost, updatedevice, devicecreate)')
            ]),
        'fini': lambda args:
            List(body=[Call("nvtxRangePop")]),
        'num-devices': lambda args, retobj:
#            Call('acc_get_num_devices', args, retobj=retobj),
            None,
        'set-device': lambda args:
            #Call('acc_set_device_num', args),
            None,
        # Pragmas
        'atomic': None, # c.Pragma('acc atomic update'),
        'map-enter-to': lambda i, j:
            None, #c.Pragma('acc enter data copyin(%s%s)' % (i, j)),
        'map-enter-to-wait': lambda i, j, k:
            None, #(c.Pragma('acc enter data copyin(%s%s) async(%s)' % (i, j, k)),
            # c.Pragma('acc wait(%s)' % k)),
        'map-enter-alloc': lambda i, j:
            None, #c.Pragma('acc enter data create(%s%s)' % (i, j)),
        'map-present': lambda i, j:
            None, #c.Pragma('acc data present(%s%s)' % (i, j)),
        'map-wait': lambda i:
            None, #c.Pragma('acc wait(%s)' % i),
        'map-update': lambda i, j:
            None, #c.Pragma('acc update self(%s%s)' % (i, j)),
        'map-update-host-if': lambda i, j, k:
            None, #c.Pragma('acc update self(%s%s) if(%s)' % (i, j, k)),
        'map-update-host': lambda i, j:
            None, #c.Pragma('acc update self(%s%s)' % (i, j)),
        'map-update-host-async': lambda i, j, k:
            None, #c.Pragma('acc update self(%s%s) async(%s)' % (i, j, k)),
        'map-update-host-async-if': lambda i, j, k, l:
            None, #c.Pragma('acc update self(%s%s) async(%s) if(%s)' % (i, j, k, l)),
        'map-update-device': lambda i, j:
            None, #c.Pragma('acc update device(%s%s)' % (i, j)),
        'map-update-device-async': lambda i, j, k:
            None, #c.Pragma('acc update device(%s%s) async(%s)' % (i, j, k)),
        'map-update-device-async-if': lambda i, j, k, l:
            None, #c.Pragma('acc update device(%s%s) async(%s) if(%s)' % (i, j, k, l)),
        'map-release': lambda i, j:
            None, #c.Pragma('acc exit data delete(%s%s)' % (i, j)),
        'map-release-if': lambda i, j, k:
            None, #c.Pragma('acc exit data delete(%s%s) if(%s)' % (i, j, k)),
        'map-exit-delete': lambda i, j:
            None, #c.Pragma('acc exit data delete(%s%s)' % (i, j)),
        'map-exit-delete-if': lambda i, j, k:
            None, #c.Pragma('acc exit data delete(%s%s) if(%s)' % (i, j, k)),
        'memcpy-to-device': lambda i, j, k:
            Call('acc_memcpy_to_device', [i, j, k]),
        'memcpy-to-device-wait': lambda i, j, k, l:
            Lambda(body=[Call('acc_memcpy_to_device_async', [i, j, k, l]),
                       Call('acc_wait', [l])]),
        'device-get':
            #Call('acc_get_device_num'),
            Call('max', (0, 0,)),
        'device-alloc': lambda i, *a, retobj=None:
            Conditional(CondEq(VOID(retobj, '*'), 0), List(body=[c.Statement(f'printf("allocating %d bytes for {retobj}\\n", {ccode(i)})'), CudaChecked(Call('cudaMalloc', (VOID(Byref(retobj), '**'), i,)))])),
        'device-free': lambda i, *a:
            #Call('acc_free', (i,))
            CudaChecked(Call('cudaFree', (i,))),
        'host-alloc': lambda i, j, k: # this isn't really 'host', it's 'high bandwidth memory'
            CudaChecked(Call("cudaMallocHost", (i, k,))),
        'host-free': lambda i:
            CudaChecked(Call("cudaFreeHost", (i,))),
        'wait-event': lambda i, j:
            CudaChecked(Call("cudaStreamWaitEvent", (i, j,))),
        'create-event': lambda i:
            CudaChecked(Call("cudaEventCreateWithFlags", (Byref(i), "cudaEventDisableTiming"))),
        'destroy-event': lambda i:
            CudaChecked(Call("cudaEventDestroy", (i,))),
        'record-event': lambda i, j:
            CudaChecked(Call("cudaEventRecord", (i, j)))
    }

    Region = OmpRegion
    HostIteration = OmpIteration  # Host parallelism still goes via OpenMP
    DeviceIteration = DeviceCudaIteration

    @classmethod
    def _map_alloc(cls, f, imask=None, condition=None):
        return CudaAlloc(f, imask, condition)

    @classmethod
    def _map_to_wait(cls, f, imask=None, qid=None):
        return None
        #return CudaTransfer(f, imask, None, CudaTransferDirection.H2D)
        #return PragmaTransfer(cls.mapper['map-enter-to-wait'], f, imask, qid)

    @classmethod
    def _map_present(cls, f, imask=None):
        return PragmaTransfer(cls.mapper['map-present'], f, imask)

    @classmethod
    def _map_wait(cls, qid=None):
        return Pragma(cls.mapper['map-wait'], qid)

    @classmethod
    def _map_update_host(cls, f, imask=None, condition=None):
        return CudaTransfer(f, imask, condition, CudaTransferDirection.D2H)

    @classmethod
    def _map_update_device(cls, f, imask=None, condition=None):
        return CudaTransfer(f, imask, condition, CudaTransferDirection.H2D)

    @classmethod
    def _map_delete(cls, f, imask=None, devicerm=None):
        return CudaDealloc(f, imask, devicerm)

    @classmethod
    def _map_update_host_async(cls, f, imask=None, qid=None, condition=None):
        return CudaTransfer(f, imask, condition, CudaTransferDirection.D2H, stream=qid)

    @classmethod
    def _map_update_device_async(cls, f, imask=None, qid=None, condition=None):
        return CudaTransfer(f, imask, condition, CudaTransferDirection.H2D, stream=qid)

    @classmethod
    def _map_wait_event(cls, e, stream=None):
        if stream is None:
            stream = e.stream if e.stream is not None else 0
        return List(body=[cls.mapper['wait-event'](stream, e.event)])

    @classmethod
    def _map_wait_recreate_event(cls, e, stream=None):
        if stream is None:
            stream = e.stream if e.stream is not None else 0
        return List(body=[cls.mapper['wait-event'](stream, e.event),
                          cls.mapper['destroy-event'](e.event),
                          cls.mapper['create-event'](e.event),
                          ])

    @classmethod
    def _map_create_event(cls, e):
        return cls.mapper['create-event'](e.event)

    @classmethod
    def _map_recreate_event(cls, e):
        return List(body=[cls.mapper['destroy-event'](e.event),
                          cls.mapper['create-event'](e.event),
                          ])

    @classmethod
    def _map_fire_event(cls, e, stream=None):
        if stream is None:
            stream = e.stream if e.stream is not None else 0
        return cls.mapper['record-event'](e.event, stream)


class DeviceCudaizer(PragmaDeviceAwareTransformer):

    lang = CudaBB
    DeviceIteration = lang.DeviceIteration

    count = 0
    def _extract_kernels(self, candidates, nthreads=None):
        assert candidates

        #root, collapsable = self._select_candidates(candidates)
        #ncollapsable = len(collapsable)

        root = candidates[0]
        if self._is_offloadable(root):
            kernel_name = "kernel%s" % (self.count)
            #body = self.DeviceIteration(gpu_fit=self.gpu_fit,
            #                            ncollapse=0,
            #                            **root.args)
            kernel, extracted_iterators = self._make_cuda_kernel(kernel_name, root)
            # find the non-derived dimensions we're iterating over, since the dimension
            # list for an Iteration includes the original dimension and the derived version
            kdims = [next(filter(lambda x: x.is_Derived == False, c.dimensions)).symbolic_size for c in extracted_iterators][:3]

            # If we don't find any non-derived dimensions, use the derived ones I guess?
            if len(kdims) == 0:
                kdims = [c.dimensions[0].symbolic_size for c in extracted_iterators][:3]

            if len(kdims) == 0:
                return root, None, None

            # The GPU wants to iterate over the whole problem space for memory alignment reasons;
            # filtering out unwanted points in the space on the GPU is very approximately zero-cost
            kgrid = kdims.copy()

            kthread = [1] * len(kdims)
            # CUDA has hard limits on the maximum size of a thread block and of grid dimensions
            # so adjust accordingly
            #if kgrid[1] > 1024:
            #    kthread[1] = 1024
            #    kgrid[1] = kgrid[1] / 1024

            #if len(kthread) > 2:
            #    kthread[0] = 1
            #    kthread[1] = 1
            #    kthread[2] = 64
            #    kgrid[0] = kgrid[0] / 1
            #    kgrid[1] = kgrid[1]# / 2# / 8
            #    kgrid[2] = kgrid[2] / 64


            partree = CudaCall(kernel_name, kgrid, kthread, kernel.parameters, stream=KernelStream())
            # Make sure that the enclosing function knows we need the full size of the Functions
            partree.expr_symbols = as_tuple(flatten((partree.expr_symbols, kdims)))
            partree = List(body=[
                #Call("setupGrid", (kthread, kgrid, kgrid[0], kgrid[1] if len(kgrid) > 1 else 1, kgrid[2] if len(kgrid) > 2 else 1)),
                partree,
            ])

            self.count = self.count + 1

            return root, partree, kernel

        elif not self.par_disabled:
            # Resort to host parallelism
            root, partree = super()._make_partree(candidates, nthreads)
            return root, partree, None

        else:
            return root, None, None

    def _make_parallel(self, iet):
        mapper = {}
        kernels = []

        # the _cudaChecked function is defined in devito_cuda.cuh

        for tree in retrieve_iteration_tree(iet, mode='superset'):
            # Get the parallelizable Iterations in `tree`
            candidates = filter_iterations(tree, key=self.key)
            if not candidates:
                continue

            # Outer parallelism
            root, partree, kernel = self._extract_kernels(candidates)
            if partree is None or root in mapper:
                continue

            mapper[root] = partree
            kernels.append(kernel)


        iet = Transformer(mapper).visit(iet)
        attrs = {'efuncs': kernels, 'includes': self.lang['headers']}

        # Also, insert CUDA events
        # sync_spots = FindNodes(SyncSpot).visit(iet)
        # if not sync_spots:
        #     return iet, attrs

        # subs = {}
        # for n in sync_spots:
        #     cuda_waits = [x for x in n.sync_ops if isinstance(x, CudaWaitEvent)]
        #     if cuda_waits:
        #         subs[n] = (n, [self.lang._map_wait_recreate_event(w) for w in cuda_waits])

        # iet = Transformer(subs).visit(iet)

        #iet = iet._rebuild(body=iet.body._rebuild(objs=iet.body.objs + (HostStream(), MemCopyStream(),)))
        return iet, attrs

    def _make_nested_partree(self, partree):
        if isinstance(partree, Callable) or isinstance(partree.root, self.DeviceIteration):
            # no-op for now
            return partree
        else:
            return super()._make_nested_partree(partree)

    def _make_cuda_kernel(self, name, body):
        # Need to extract the outermost Iterations, as they'll be handled by the GPU hardware
        dim_iter = []
        #node = body#self._make_reductions(body)
        #while isinstance(node, Iteration) and len(dim_iter) < 3:
        #    iter = node
        #    dim_iter.append(iter)
        #    if len(node.nodes) == 1:
        #        node = node.nodes[0]
        #        iter.nodes = []
        #    else:
        #        node = node.nodes
        #        break

        # Find the iterators we consider eligible for being the GPU grid dimensions
        iterations = FindNodes(Iteration).visit(body)
        possible_iter_dimensions = list(OrderedDict.fromkeys([x.dim for x in iterations if not x.dim.name.startswith("par_dim")]))
        grouped_iters = [(x, list(OrderedDict.fromkeys([i for i in iterations if i.dim == x]))) for x in possible_iter_dimensions]
        valid_dims = list(filter(lambda i: len(set([z.limits for z in i[1]])) == 1, grouped_iters))
        if len(valid_dims) > 3:
            valid_dims = valid_dims[0:3]

        iet = IterationExtractor([d[0] for d in valid_dims]).visit(body)
        
        # Now, generate the iteration dimension variables from the blockIdx/threadIdx
        dim_vars = ["x", "y", "z"]
        kernel = []
        args = set()
        for v in range(0, len(valid_dims)):
            dim, iters = valid_dims[v]
            limits = iters[0].limits
            symbols = flatten([i.expr_symbols for i in iters])
            # FIXME: this should be a declaration, not an assignment - figure out how to get the iteration variables
            # out of the signature?
            kernel.append(c.Initializer(c.Value('int', dim.name), "blockDim.%s * blockIdx.%s + threadIdx.%s" % (dim_vars[v], dim_vars[v], dim_vars[v])))
            args = args.union(symbols)
            # Add the iteration conditions
            kernel.append(c.If("%s < %s || %s >= %s" % (dim.name, str(limits[0]), dim.name, str(limits[1])), c.Statement("return")))

        # Add the iteration body
        kernel.extend(as_tuple(iet))

        # Remove the original iteration variables from the signature

        cuda_callable = CudaCallable(name=name, body=kernel, parameters=args, defines=[x[0] for x in valid_dims])

        return (cuda_callable, list([x[1][0] for x in valid_dims]))

    def _make_reductions(self, partree):
        if not any(i.is_ParallelAtomic for i in partree.collapsed):
            return partree

        exprs = [i for i in FindNodes(Expression).visit(partree) if i.is_reduction]
        reductions = [(i.output, i.operation) for i in exprs]

        test0 = all(not i.is_Indexed for i, _ in reductions)

        if test0:
            # Implement reduction
            mapper = {partree.root: partree.root._rebuild(reduction=reductions)}
        elif all(i is OpInc for _, i in reductions):
            # Use atomic increments
            mapper = {i: i._rebuild(pragmas=self.lang['atomic']) for i in exprs}
        else:
            raise NotImplementedError

        partree = Transformer(mapper).visit(partree)

        return partree

class IterationExtractor(Visitor):
    def __init__(self, dims):
        super(Visitor, self).__init__()
        self._dims = dims        

    def visit_object(self, o, **kwargs):
        return o

    def visit_tuple(self, o, **kwargs):
        visited = tuple(self._visit(i, **kwargs) for i in o)
        return tuple(i for i in visited if i is not None)

    visit_list = visit_tuple

    def visit_Iteration(self, o, **kwargs):
        if o.dim in self._dims:
            return List(body=self._visit(o.nodes, **kwargs))
        
        else:
            children = [self._visit(i, **kwargs) for i in o.children]
            return o._rebuild(*children, **o.args_frozen)
    
    def visit_Node(self, o, **kwargs):
        children = [self._visit(i, **kwargs) for i in o.children]
        return o._rebuild(*children, **o.args_frozen)

    
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
        self.gpu_fit = options['gpu-fit']

    def _alloc_local_array_on_high_bw_mem(self, site, obj, storage, devicerm=None):
        """
        Allocate a local Array in the device high bandwidth memory.
        """
        # Create a local static pointer so we can persist between runs
        decl = Definition(obj, initvalue="nullptr", prefix="static")
        doalloc = self.lang['device-alloc']
        dofree = self.lang['device-free']

        nbytes = SizeOf(obj._C_typedata)*obj.size
        init = doalloc(nbytes, None, retobj=obj._C_symbol)
        allocs = (init, ) if isinstance(init, Call) and init.retobj == obj else (decl, init)

        free = dofree(obj._C_name, None)

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
            definition = (decl)

        frees = obj._C_free

        storage.update(obj, site, objs=definition, frees=frees)

    def _alloc_array_on_low_lat_mem(self, site, obj, storage):
        """
        Allocate an Array in the low latency memory.
        """
        shape = "".join("[%s]" % ccode(i) for i in obj.symbolic_shape)
        alignment = self.lang['aligned'](obj._data_alignment)
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

        memptr = VOID(Byref(obj._C_symbol), '**')
        alignment = obj._data_alignment
        nbytes = SizeOf(obj._C_typedata)*obj.size
        alloc = self.lang['host-alloc'](memptr, alignment, nbytes)

        free = self.lang['host-free'](obj._C_symbol)

        storage.update(obj, site, allocs=(decl, alloc), frees=free)

    def _alloc_mapped_array_on_high_bw_mem(self, site, obj, storage, *args):
        """
        Allocate a mapped Array in the host high bandwidth memory.
        """
        static_decl = Definition(obj, initvalue="nullptr", prefix="static")
        decl = Definition(obj, initvalue="nullptr")
        # Allocating a mapped Array on the high bandwidth memory requires
        # multiple statements, hence we implement it as a generic Callable
        # to minimize code size, since different arrays will ultimately be
        # able to reuse the same abstract Callable

        memptr = VOID(Byref(obj._C_symbol), '**')
        alignment = obj._data_alignment
        nbytes = SizeOf(obj._C_typedata)
        alloc0 = self.lang['host-alloc'](memptr, alignment, nbytes)

        nbytes_param = Symbol(name='nbytes', dtype=np.uint64, is_const=True)
        nbytes_arg = SizeOf(obj.indexed._C_typedata)*obj.size

        ffp1 = FieldFromPointer(obj._C_field_data, obj._C_symbol)
        memptr = VOID(Byref(ffp1), '**')
        alloc1 = self.lang['host-alloc'](memptr, alignment, nbytes_param)

        ffp2 = FieldFromPointer(obj._C_field_device_data, obj._C_symbol)
        memptr2 = VOID(Byref(ffp2), '*')
        alloc2 = self.lang['device-alloc'](nbytes_param, retobj=ffp2)

        ffp0 = FieldFromPointer(obj._C_field_nbytes, obj._C_symbol)
        init0 = DummyExpr(ffp0, nbytes_param)
        init1 = DummyExpr(ffp1, 0)
        init2 = DummyExpr(ffp2, 0)

        free0 = self.lang['host-free'](ffp1)

        free1 = self.lang['device-free'](ffp2)

        free2 = self.lang['host-free'](obj._C_symbol)

        ret = Return(obj._C_symbol)

        alloc_name = self.sregistry.make_name(prefix='alloc')
        body = (decl, alloc0, init0, init1, init2, alloc1, alloc2, ret)
        #body = (decl, alloc0, alloc1, init, ret)
        efunc0 = make_callable(alloc_name, body, retval=obj._C_typename)
        assert len(efunc0.parameters) == 1  # `nbytes_param`

        free_name = self.sregistry.make_name(prefix='free')
        efunc1 = make_callable(free_name, (free0, free1, free2))
        #efunc1 = make_callable(name, (free0, free2))
        assert len(efunc1.parameters) == 1  # `obj`
        alloc = List(body=[static_decl, Conditional(CondOr(CondEq(VOID(obj._C_symbol, '*'), 0), CondNe(nbytes_arg, ffp0)),
                            List(body=[
                                    Conditional(CondNe(VOID(obj._C_symbol, '*'), 0), Block(body=[c.Statement(f'printf("resizing {obj._C_symbol} from %d to %d bytes\\n", {ccode(ffp0)}, {ccode(nbytes_arg)})'),
                                                                                                 Call(free_name, obj),
                                                                                                 c.Assign(obj._C_symbol, 0)])),
                                    Call(alloc_name, nbytes_arg, retobj=obj, declares=False)]))])
        free = List(body=[Conditional(DeviceRM(), List(body=[c.Statement(f'printf("deleting local storage for {obj._C_symbol}\\n")'), Call(free_name, obj)]))])

        storage.update(obj, site, allocs=alloc, frees=free, efuncs=(efunc0, efunc1))

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

        memptr = VOID(Byref(obj._C_symbol), '**')
        alignment = obj._data_alignment
        nbytes = SizeOf(Keyword('%s*' % obj._C_typedata))*obj.dim.symbolic_size
        alloc0 = self.lang['host-alloc'](memptr, alignment, nbytes)

        free0 = self.lang['host-free'](obj._C_symbol)

        # The pointee Array
        pobj = IndexedPointer(obj._C_symbol, obj.dim)
        memptr = VOID(Byref(pobj), '**')
        nbytes = SizeOf(obj._C_typedata)*obj.array.size
        alloc1 = self.lang['host-alloc'](memptr, alignment, nbytes)

        free1 = self.lang['host-free'](pobj)

        # Dump
        if obj.dim is self.sregistry.threadid:
            storage.update(obj, site, allocs=(decl, alloc0), frees=free0,
                           pallocs=(obj.dim, alloc1), pfrees=(obj.dim, free1))
        else:
            storage.update(obj, site, allocs=(decl, alloc0, alloc1), frees=(free0, free1))

    def _map_array_on_high_bw_mem(self, _site, _obj, _storage):
        """
        Map an Array already defined in the host memory in to the device high
        bandwidth memory.
        """

        # When using CUDA we allocate everything in a device-visible manner
        return

    def _map_function_on_high_bw_mem(self, site, obj, storage, devicerm, read_only=False, devicecreate=None, updatehost=None, updatedevice=None):
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
            mmap = [#self.lang._map_alloc(obj, condition=devicecreate),
                    self.lang._map_update_device(obj, condition=CondOr(devicecreate, updatedevice))]
        else:
            mmap = self.lang._map_to(obj)

        if read_only is False:
            unmap = [self.lang._map_update_host(obj, condition=CondOr(updatehost, devicerm)),
                     self.lang._map_release(obj, devicerm=devicerm),
            ]
        else:
            unmap = self.lang._map_delete(obj, devicerm=devicerm)

        storage.update(obj, site, maps=mmap, unmaps=unmap)

    def _dump_transfers(self, iet, storage):
        mapper = {}
        for k, v in storage.items():
            if v.maps or v.unmaps:
                mapper[iet.body] = iet.body._rebuild(maps=flatten(v.maps),
                                                     unmaps=flatten(v.unmaps))

        processed = Transformer(mapper, nested=True).visit(iet)

        return processed

    @iet_visit
    def derive_transfers(self, iet):
        """
        Collect all symbols that cause host-device data transfer, distinguishing
        between reads and writes.
        """

        def needs_transfer(f):
            return (f._mem_mapped and
                    not isinstance(f, AliasFunction) and
                    is_on_device(f, self.gpu_fit))

        writes = set()
        reads = set()
        for i, v in MapExprStmts().visit(iet).items():
            if not any(isinstance(j, self.lang.DeviceIteration) for j in v) and \
               not isinstance(i, DeviceCall) and \
               not isinstance(iet, DeviceFunction):
                # Not an offloaded Iteration tree
                continue

            writes.update({w for w in i.writes if needs_transfer(w)})
            reads.update({f for f in i.functions
                          if needs_transfer(f) and f not in writes})

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
                    self._map_function_on_high_bw_mem(iet, i, storage, devicerm, devicecreate = devicecreate, updatehost = updatehost, updatedevice = updatedevice)
            for i in filter_sorted(reads - writes):
                if i.is_Array:
                    self._map_array_on_high_bw_mem(iet, i, storage)
                else:
                    self._map_function_on_high_bw_mem(iet, i, storage, devicerm, True, devicecreate = devicecreate, updatehost = updatehost, updatedevice = updatedevice)

            iet = self._dump_transfers(iet, storage)

            return iet, {}

        return _place_transfers(iet, mapper=kwargs['mapper'])

    @iet_pass
    def place_cuda_casts(self, iet, **kwargs):
        # Don't generate unnecessary casts in the entry function and in CUDA kernels
        if not isinstance(iet, CudaCallable) and not isinstance(iet, EntryFunction):
            # Candidates
            indexeds = FindSymbols('indexeds|indexedbases').visit(iet)

            # Create Function -> n-dimensional array casts
            # E.g. `float (*u)[.] = (float (*)[.]) u_vec->data`
            # NOTE: a cast is needed only if the underlying data object isn't already
            # defined inside the kernel, which happens, for example, when:
            # (i) Dereferencing a PointerArray, e.g., `float (*r0)[.] = (float(*)[.]) pr0[.]`
            # (ii) Declaring a raw pointer, e.g., `float * r0 = NULL; *malloc(&(r0), ...)
            defines = set(FindSymbols('defines').visit(iet))
            bases = sorted({i.base for i in indexeds}, key=lambda i: i.name)
            casts = [self.lang.PointerCast(i.function, obj=i) for i in bases
                    if i not in defines]

            # Incorporate the newly created casts
            if casts:
                iet = iet._rebuild(body=iet.body._rebuild(casts=casts))

        return iet, {}

    @iet_pass
    def tidy_up(self, iet, **kwargs):
        class Tidier(Visitor):
            """
            A basic Visitor that performs rudimentary tidying up.

            Currently limited to coalescing consecutive identical conditionals
            """
            def visit_object(self, o):
                return o
            def visit_list(self, o):
                return o
            def visit_tuple(self, o):
                return o
            def visit_Collection(self, o):
                return o
            def visit_List(self, o):
                return o

        return iet, {}#Tidier().visit(iet), {}

    def process(self, graph):
        """
        Apply the `place_transfers`, `place_definitions` and `place_casts` passes.
        """
        mapper = self.derive_transfers(graph)
        self.place_transfers(graph, mapper=mapper)
        self.place_definitions(graph)
        self.place_cuda_casts(graph)
        self.tidy_up(graph)


class CudaHostFuncCall(AsyncCall):
    def __init__(self, name, arguments=None, retobj=None, is_indirect=False,
                 cast=False, writes=None, types=None, declares=True, stream=None):
        super().__init__(name, arguments=arguments, retobj=retobj, is_indirect=is_indirect,
                         cast=cast, writes=writes, types=types, declares=declares)
        self.stream = stream

class CudaHostFuncCallable(AsyncCallable):
    pass

class KernelStream(CudaStream):
    def __init__(cls):
        super().__init__("kernel_stream")

    def __new__(cls):
        return super().__new__(cls, name="kernel_stream")
class MemCopyStream(CudaStream):
    def __init__(cls):
        super().__init__("memcpy_stream")

    def __new__(cls):
        return super().__new__(cls, name="memcpy_stream")

class HostStream(CudaStream):
    def __init__(cls):
        super().__init__("host_stream")
    def __new__(cls):
        return super().__new__(cls, name="host_stream")

class CudaOrchestrator(Orchestrator):
    lang = CudaBB

    _memcpy_stream = MemCopyStream()
    _host_stream = HostStream()
    _kernel_stream = KernelStream()

    def _make_waitlock(self, iet, sync_ops):
        waitloop = List(
            header=c.Comment("Wait for `%s` to be copied to the host" %
                             ",".join(s.function.name for s in sync_ops)),
            body=[self.lang._map_wait_recreate_event(s, stream=self._kernel_stream) for s in sync_ops],
            footer=c.Line()
        )

        #trigger = List(
        #    header = c.Comment("Let the background stream know we're done with `%s`" %
        #                        ",".join(s.function.name for s in sync_ops)),
        #    body=[self.lang._map_fire_event(s, stream=self._kernel_stream) for s in sync_ops],
        #    footer=c.Line()
        #)
        #iet = List(body=flatten([(waitloop,), iet.body, trigger]))
        iet = List(body=flatten([(waitloop,), iet.body]))
        return iet, []

    def _make_waitevent(self, iet, sync_ops):
        waitloop = List(
            header=c.Comment("Wait for `%s` to be copied to the host" %
                             ",".join(s.function.name for s in sync_ops)),
            body=[self.lang._map_wait_recreate_event(s, stream=self._kernel_stream) for s in sync_ops],
            footer=c.Line()
        )

        iet = List(body=(waitloop,) + iet.body)

        return iet, []

    def _make_releaselock(self, iet, sync_ops):
        preactions = []
        preactions.extend(self.lang._map_fire_event(s) for s in sync_ops)

        iet = List(
            header=c.Comment("Release lock(s) as soon as possible"),
            body=preactions + [iet]
        )

        return iet, []

    def _make_withlock(self, iet, sync_ops):
        qid = QueueID()

        preactions = [c.Comment("Block the copy until it's safe"),
                      BlankLine]

        # these should run on the memcpy stream, so it needs to wait
        preactions.extend([self.lang._map_wait_event(s, stream=self._memcpy_stream) for s in sync_ops])
        # and the main kernel stream should mark this as the appropriate place for it to start
        preactions.extend([self.lang._map_fire_event(s, stream=self._kernel_stream) for s in sync_ops])
        # then recreate the event
        preactions.extend([self.lang._map_recreate_event(s) for s in sync_ops])
        preactions.extend([self.lang._map_update_host_async(s.function, qid=self._memcpy_stream) for s in sync_ops])
        preactions.extend([self.lang._map_fire_event(s, stream=self._memcpy_stream) for s in sync_ops])

        # these should run on the host stream (not default stream)
        preactions.extend([self.lang._map_wait_recreate_event(s, stream=self._host_stream) for s in sync_ops])
        postactions = [BlankLine, c.Comment("Raise the event")]
        postactions.extend([self.lang._map_fire_event(s,stream=self._host_stream) for s in sync_ops])

        # Turn `iet` into an AsyncCallable so that subsequent passes know
        # that we're happy for this Callable to be executed asynchronously
        name = self.sregistry.make_name(prefix='copy_device_to_host')
        async_body = List(body=iet.body)
        parameters = _cuda_derive_parameters(async_body)
        async_body = async_body._rebuild()
        efunc = CudaHostFuncCallable(name, async_body, parameters=parameters)

        # The corresponding AsyncCall
        body = preactions + [CudaHostFuncCall(name, efunc.parameters, stream=self._host_stream)] + postactions

        iet = List(body=body)

        return iet, [efunc]

    def _make_fetchupdate(self, iet, sync_ops):
        postactions = [self.lang._map_update_device(s.target, s.imask)
                       for s in sync_ops]

        # Turn init IET into a Callable
        name = self.sregistry.make_name(prefix='init_device')
        body = List(body=iet.body + tuple(postactions))
        parameters = derive_parameters(body)
        efunc = Callable(name, body, 'void', parameters, 'static')

        # Perform initial fetch by the main thread
        iet = List(
            header=c.Comment("Initialize data stream"),
            body=Call(name, parameters)
        )

        return iet, [efunc]

    def _make_prefetchupdate(self, iet, sync_ops):
        qid = QueueID()
        preactions = []
        preactions.extend([self.lang._map_wait_event(s, stream=self._host_stream) for s in sync_ops])        
        preactions.extend([self.lang._map_fire_event(s, stream=self._kernel_stream) for s in sync_ops])
        preactions.extend([self.lang._map_recreate_event(s) for s in sync_ops])


        
        postactions = []
        postactions.extend([self.lang._map_fire_event(s, stream=self._host_stream) for s in sync_ops])
        postactions.extend([self.lang._map_wait_recreate_event(s, stream=self._memcpy_stream) for s in sync_ops])
        postactions.extend([self.lang._map_update_device_async(s.target, qid=self._memcpy_stream) for s in sync_ops])
        postactions.extend([self.lang._map_fire_event(s, stream=self._memcpy_stream) for s in sync_ops])

        # Turn `iet` into an AsyncCallable so that subsequent passes know
        # that we're happy for this Callable to be executed asynchronously
        name = self.sregistry.make_name(prefix='prefetch_host_to_device')
        body = iet.body 
        #+ (BlankLine,) + tuple(postactions))
        parameters = _cuda_derive_parameters(body)
        efunc = CudaHostFuncCallable(name, body, parameters=parameters)

        # The corresponding AsyncCall
        iet = List(body=preactions + [CudaHostFuncCall(name, efunc.parameters, stream=self._host_stream)] + postactions)

        return iet, [efunc]

    @iet_pass
    def process(self, iet):
        sync_spots = FindNodes(SyncSpot).visit(iet)

        if not sync_spots:
            if isinstance(iet, EntryFunction):
                iet = iet._rebuild(body=List(body=[iet.body, CudaChecked(Call("cudaDeviceSynchronize", None))]))
            return iet, {}

        callbacks = OrderedDict([
            (CudaWaitEvent, self._make_waitlock),
            (CudaWithEvent, self._make_withlock),
            (CudaFireEvent, self._make_withlock),
            (CudaFetchUpdate, self._make_fetchupdate),
            (CudaPrefetchUpdate, self._make_prefetchupdate),
        ])

        # The SyncOps are to be processed in a given order
        key = lambda s: list(callbacks).index(s)

        efuncs = []
        subs = {}
        events = []
        for n in sync_spots:
            mapper = as_mapper(n.sync_ops, lambda i: type(i))
            for t in sorted(mapper, key=key):
                events.extend([e.event for e in mapper[t] if e.event is not None])
                subs[n], v = callbacks[t](subs.get(n, n), mapper[t])
                efuncs.extend(v)

        iet = Transformer(subs).visit(iet)

        events = [List(body=[Definition(e, None, None, NullPointer()),
                             self.lang.mapper['create-event'](e._C_symbol),
                            # self.lang.mapper['record-event'](e._C_symbol, KernelStream()._C_symbol)
                             ]) for e in filter_ordered(events)]
        iet = iet._rebuild(body = List(body=[events, iet.body, CudaChecked(Call("cudaDeviceSynchronize", None))]))

        return iet, {'efuncs': efuncs}

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

def _cuda_derive_parameters(iet):
    indexeds = FindSymbols('indexeds|indexedbases').visit(iet)
    bases = sorted({i.base for i in indexeds}, key=lambda i: i.name)

    return filter_ordered(flatten([derive_parameters(iet), FindSymbols('basics').visit(bases)]))

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
    extra_parameters = tuple([i for i in FindSymbols('basics').visit(casts) if i not in iet.parameters])
    track[iet.name].extra_args = extra_parameters

    fields = iet.parameters + extra_parameters
    defines = FindSymbols('defines').visit(root.body)
    ncfields, cfields = split(fields, lambda i: i in defines)

    # SharedData -- that is the data structure that will be used by the
    # main thread to pass information down to the child thread(s)
    sdata = track[iet.name].sdata = CudaSharedData(name='sdata',
                                               npthreads=1,
                                               cfields=cfields,
                                               ncfields=ncfields,
                                               pname='tsdata%d' % n)
    sbase = sdata.symbolic_base

    # Prepend the SharedData fields available upon thread activation
    #preactions = [DummyExpr(i, FieldFromPointer(i.name, sbase)) for i in ncfields]
    #preactions.append(BlankLine)
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
        #threads = b.threads._rebuild(name=name)

        # Call to `sdata` initialization Callable
        sbase = sdata.symbolic_base
        d = 0#threads.index
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
        #assert len(efuncs[n.name].parameters) == len(arguments)
        call0 = Call(efuncs[n.name].name, arguments)

        initialization.append(List(
            body=[Definition(sdata, None, None, NullPointer())],#c.Initializer(c.Value("%s*" % sdata._C_typedata, sdata._C_symbol), NullPointer())],
        ))

        # Activation
        #if threads.size == 1:
        d = 0

        activation = [c.Comment("Allocate a new block of data for this invocation"),
                      Call("posix_memalign", (VOID(Byref(sdata), '**'), 64, SizeOf(sdata._C_typedata))),
                      call0]
        activation.extend([DummyExpr(FieldFromComposite(i.name, sdata[d]), i)
                           for i in sdata.ncfields])

        activation = activation + [

           c.Statement("cudaLaunchHostFunc(%s, (cudaHostFn_t)%s, %s)" % (n.stream if n.stream is not None else 0, n.name, ccode(sbase + d))),
            ]
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
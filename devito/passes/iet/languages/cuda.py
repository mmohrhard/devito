from collections import OrderedDict
from ctypes import POINTER, c_void_p
from enum import Enum
from functools import cached_property, singledispatch
import cgen as c
from sympy import Or, simplify, Number
from devito.ir.iet.visitors import Visitor

from devito.symbolics.extended_sympy import FieldFromComposite, Null, ReservedWord
from devito.types.misc import Global, Pointer
from devito.tools.data_structures import Bunch, DefaultOrderedDict
from devito.ir.support.syncs import FetchUpdate, PrefetchUpdate, ReleaseLock, WaitLock, WithLock
from devito.ir.iet.efunc import AsyncCall, AsyncCallable, ThreadCallable
from devito.ir.iet.cuda import CudaTransferDirection, TemplateParameter
from devito.passes.iet.definitions import DeviceAwareDataManager, Storage
from devito.passes.iet.langbase import make_sections_from_imask
from devito.symbolics.printer import ccode
from devito.ir.iet.nodes import BlankLine, BusyWait, Conditional, Dereference, PointerCast, Return, Section, SyncSpot, Transfer, While
from devito.tools.utils import as_list, as_mapper, as_tuple, dtype_to_ctype, filter_sorted, flatten, split
from devito.types.parallel import CudaStream, DeviceRM, QueueID, SharedData, ThreadArray, Lock, CudaEvent
import numpy as np
from sympy import Max
from devito.logger import info, warning, debug
from devito.arch import CUDA, NVIDIAX
from devito.ir import (Call, Callable, CudaCall, CudaCallableBody, DeviceCall, DummyExpr, DPtr, EntryFunction, List, CudaCallable,
                       Block, ParallelIteration, ParallelTree, Pragma, Definition, Iteration, Node,
                       FindNodes, FindSymbols, Uxreplace, Transformer, Lambda, AddressOf, CLiteral,
                       MapExprStmts, DeviceFunction, make_callable, Expression, OpInc, CudaKernelPointerCast,
                       )
from devito.passes.iet.engine import iet_pass, iet_visit
from devito.passes.iet.definitions import DataManager
from devito.passes.iet.orchestration import Orchestrator
from devito.passes.iet.parpragma import (PragmaDeviceAwareTransformer, PragmaLangBB,
                                         PragmaTransfer, PragmaDeviceAwareDataManager)
from devito.symbolics import (Byref, DefFunction, FieldFromPointer, IndexedPointer,
                              ListInitializer, SizeOf, VOID, INT, Keyword, ccode,
                               CondEq, CondNe, CondOr)
from devito.passes.iet.languages.C import CBB
from devito.passes.iet.languages.openmp import OmpRegion, OmpIteration
from devito.passes.iet.languages.utils import make_clause_reduction
from devito.passes.iet.misc import is_on_device
from devito.ir.iet.utils import derive_parameters, retrieve_iteration_tree, filter_iterations
from devito.symbolics import Macro, cast_mapper, uxreplace
from devito.symbolics import pow_to_mul
from devito.tools import filter_ordered
from devito.types import DevicePointer, Symbol, Constant, DeviceRM, DeviceCreate, UpdateDevice, UpdateHost, Eq
from devito.types.basic import IndexedBase
from devito.types.dense import AliasFunction

__all__ = ['DeviceCudaizer', 'DeviceCudaDataManager', 'CudaOrchestrator', 'cuda_memcpy', 'cuda_eventify', 'KernelStream', 'NcclStream', 'CudaChecked']


class KernelStream(CudaStream, Global):
    def __init__(cls, *args, **kwargs):
        super().__init__("kernel_stream")

    def __new__(cls, *args):
        return super().__new__(cls, "kernel_stream")

    def __new__(cls, *args, **kwargs):
        return super().__new__(cls, "kernel_stream")

class MemCopyStream(CudaStream, Global):
    def __init__(cls, *args, **kwargs):
        super().__init__("memcpy_stream")

    def __new__(cls, *args):
        return super().__new__(cls, "memcpy_stream")

    def __new__(cls, *args, **kwargs):
        return super().__new__(cls, "memcpy_stream")

class HostStream(CudaStream, Global):
    def __init__(cls, *args, **kwargs):
        super().__init__("host_stream")

    def __new__(cls, *args):
        return super().__new__(cls, "host_stream")

    def __new__(cls, *args, **kwargs):
        return super().__new__(cls, "host_stream")

class NcclStream(CudaStream, Global):
    def __init__(cls, *args, **kwargs):
        super().__init__("nccl_stream")

    def __new__(cls, *args):
        return super().__new__(cls, "nccl_stream")

    def __new__(cls, *args, **kwargs):
        return super().__new__(cls, "nccl_stream")

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
    def name(self):
        return self.function._C_name
    
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
    def __init__(self, arguments=None, name=None):
        super().__init__("CudaChecked", arguments=tuple(flatten([arguments])))

    @cached_property
    def expr_symbols(self):
        return flatten([x.expr_symbols for x in flatten(self.arguments)])

class NullPointer(ReservedWord):
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls, "nullptr")

class JitifyCache(Global):
    @property
    def _C_typename(self):
        return "jitify::JitCache"    
    def __init__(cls, name, *args, **kwargs):
        super().__init__(name, dtype=c_void_p)

    def __new__(cls, *args):
        return super().__new__(cls, "kernel_cache")

    def __new__(cls, *args, **kwargs):
        return super().__new__(cls, "kernel_cache")
    
class JitifyProgram(Global):
    @property
    def _C_typename(self):
        return "jitify::Program"    
    def __init__(cls, name, *args, **kwargs):
        super().__init__(name, dtype=c_void_p)

    def __new__(cls, name, *args):
        return super().__new__(cls, name)

    def __new__(cls, name, *args, **kwargs):
        return super().__new__(cls, name)
    
class CudaBB(PragmaLangBB):

    mapper = {
        # Misc
        'name': 'CUDA',
        'headers': ['cuda.h', 'cuda_runtime_api.h', 'nvtx3/nvToolsExt.h', 'stdio.h', 'assert.h', 'devito/devito_cuda.cuh', 'nccl.h', 'devito/jitify.hpp'],
        'global-decls': [
            Definition(HostStream(), initvalue = "nullptr", prefix="static"),
            Definition(MemCopyStream(), initvalue = "nullptr", prefix="static"),
            Definition(KernelStream(), initvalue = "nullptr", prefix="static"),
            Definition(NcclStream(), initvalue = "nullptr", prefix="static")
        ],
        # Platform mapping
        CUDA: None,
        NVIDIAX: None,
        # Runtime library
        'aligned': lambda i:
            '__attribute__((aligned(%d)))' % i,
        'init': lambda args:
            List(body=[Conditional(CondEq(HostStream(), NullPointer()), Call("cudaStreamCreateWithFlags", (Byref(HostStream()), "cudaStreamNonBlocking"))),
                       Conditional(CondEq(MemCopyStream(), NullPointer()), Call("cudaStreamCreateWithFlags", (Byref(MemCopyStream()), "cudaStreamNonBlocking"))),
                       Conditional(CondEq(KernelStream(), NullPointer()), Call("cudaStreamCreateWithFlags", (Byref(KernelStream()), "cudaStreamNonBlocking"))),
                       Conditional(CondEq(NcclStream(), NullPointer()), Call("cudaStreamCreateWithFlags", (Byref(NcclStream()), "cudaStreamNonBlocking"))),

                       Definition(JitifyCache("kernel_cache"), prefix="static"),
                       Definition(JitifyProgram("program"), initvalue=Call("kernel_cache.program", ("_cudaKernels", 0))),
                       Call("nvtxRangePush", ("__FUNCTION__", )),
            ]),
        'fini': lambda args:
            List(body=[Call("nvtxRangePop")]),
        'num-devices': lambda args, retobj:
            Block(body=[c.Initializer(c.Value("int", "_num_devices"), 0),
                        Call("cudaGetDeviceCount", (INT(Byref(retobj), '*')))]),
        'set-device': lambda device:
            CudaChecked(Call("cudaSetDevice", (device,))),
        # Pragmas
        'atomic': None, # CUDA doesn't use a pragma for this
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
            # calls a helper function since we expect a return value
            Call('_cudaGetCurrentDevice'),
        'device-alloc': lambda i, *a, retobj=None:
            CudaChecked(Call('cudaMalloc', (VOID(Byref(retobj), '**'), i,))),
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
    def _map_release(cls, f, imask=None, devicerm=None):
        return CudaDealloc(f, imask, devicerm)
    
    @classmethod
    def _map_update_host_async(cls, f, imask=None, qid=None, condition=None):
        return CudaTransfer(f, imask, condition, CudaTransferDirection.D2H, stream=qid)

    @classmethod
    def _map_update_device_async(cls, f, imask=None, qid=None, condition=None):
        return CudaTransfer(f, imask, condition, CudaTransferDirection.H2D, stream=qid)

    @classmethod
    def _map_wait_event(cls, e, stream=None):
        stream = stream if stream is not None else 0
        return List(body=[cls.mapper['wait-event'](stream, e.handle)])

    @classmethod
    def _map_wait_recreate_event(cls, e, stream=None):
        stream = stream if stream is not None else 0
        return List(body=[cls.mapper['wait-event'](stream, e.handle),
                          cls.mapper['destroy-event'](e.handle),
                          cls.mapper['create-event'](e.handle),
                          ])

    @classmethod
    def _map_create_event(cls, e):
        return cls.mapper['create-event'](e.handle)

    @classmethod
    def _map_recreate_event(cls, e):
        return List(body=[cls.mapper['destroy-event'](e.handle),
                          cls.mapper['create-event'](e.handle),
                          ])

    @classmethod
    def _map_fire_event(cls, e, stream=None):
        stream = stream if stream is not None else 0
        return cls.mapper['record-event'](e.handle, stream)

    @classmethod
    def _get_num_devices(cls, platform):
        ngpus = Symbol(name='_num_gpus')
        return ngpus, List(body=[
            c.Initializer(c.Value("int", "_num_gpus"), 0),
            Call("cudaGetDeviceCount", (INT(Byref(ngpus), '*')))
        ])

class CudaAtomicExpression(Expression):
    def __init__(self, expr, pragmas=None, init=None, operation=None):
        super().__init__(expr, pragmas, init, operation, True)

class DeviceCudaizer(PragmaDeviceAwareTransformer):

    lang = CudaBB
    DeviceIteration = lang.DeviceIteration

    count = 0
    def _extract_kernels(self, candidates, nthreads=None):
        assert candidates

        root = candidates[0]
        if self._is_offloadable(root):
            kernel_name = f"{self.kernel_basename}{self.count}"

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

            # the grid/threads are (for now) set up in some C++ code from a header
            kthread = [1] * len(kdims)

            partree = CudaCall(kernel_name, kgrid, kthread, preferred_block=kernel.preferred_block, preferred_sub_block=kernel.preferred_sub_block, arguments=kernel.parameters, kernel=kernel, stream=KernelStream())
            # Make sure that the enclosing function knows we need the full size of the Functions
            partree.expr_symbols = as_tuple(flatten((partree.expr_symbols, kdims)))

            self.count = self.count + 1

            return root, partree, [kernel]

        elif not self.par_disabled:
            # Resort to host parallelism
            root, partree = super()._make_partree(candidates, nthreads)
            return root, partree, None

        else:
            return root, None, None

    def _make_parallel(self, iet):
        mapper = {}
        kernels = []

        # Name kernels according to the name of the EntryFunction by default so that
        # profiling multiple operators in a single Nsight run produces more
        # meaningful summary data
        try:
            self.kernel_basename = FindNodes(EntryFunction).visit(iet)[0].name + "_kernel"
        except IndexError:
            self.kernel_basename = "kernel"

        for tree in retrieve_iteration_tree(iet, mode='superset'):
            # Get the parallelizable Iterations in `tree`
            candidates = filter_iterations(tree, key=self.key)
            if not candidates:
                continue

            # Outer parallelism
            root, partree, kernels_gen = self._extract_kernels(candidates)
            if partree is None or root in mapper:
                continue

            mapper[root] = partree
            if kernels_gen:
                kernels.extend(kernels_gen)


        iet = Transformer(mapper).visit(iet)
        attrs = {'efuncs': kernels, 'includes': self.lang['headers'], 'globals': self.lang['global-decls']}

        return iet, attrs

    def _make_nested_partree(self, partree):
        if isinstance(partree, Callable) or isinstance(partree.root, self.DeviceIteration):
            # no-op for now
            return partree
        else:
            return super()._make_nested_partree(partree)

    def _make_cuda_kernel(self, name, body):
        body = realign_iet(body)
        # Find the iterators we consider eligible for being the GPU grid dimensions
        iterations = list([i for i in FindNodes(Iteration).visit(body) if i.is_ParallelRelaxed])
        possible_iter_dimensions = list(OrderedDict.fromkeys([x.dim for x in iterations if not x.dim.name.startswith("par_dim")]))
        grouped_iters = [(x, list(OrderedDict.fromkeys([i for i in iterations if i.dim == x]))) for x in possible_iter_dimensions]
        valid_dims = list(filter(lambda i: len(set([z.limits for z in i[1]])) == 1, grouped_iters))
        if len(valid_dims) > 3:
            valid_dims = valid_dims[0:3]

        iet = IterationExtractor([d[0] for d in valid_dims]).visit(body)

        # if any of the iterations we're extracting have the atomic flag set, we need to force all reductions
        # in the kernel to be atomic regardless of any inner iterations
        force_atomic = any(i.is_ParallelAtomic for i in flatten([d[1] for d in valid_dims]))
        iet = self._make_reductions(iet, force_atomic=force_atomic)
        # replace any atomic ops
        exprs = [e for e in FindNodes(Expression).visit(iet) if e.is_atomic]
        mapper = dict([(i, CudaAtomicExpression(i.expr, i.pragmas, i.init, i.operation),) for i in exprs])
        iet = Transformer(mapper).visit(iet)

        # Now, generate the iteration dimension variables from the blockIdx/threadIdx
        # These end up reversed because warp thread order in CUDA for >1D is column-major
        # and we want contiguous memory access
        dim_vars = ["x", "y", "z"]

        kernel = []
        args = set()

        sub_blocks = [2] * (len(valid_dims) - 1)

        sub_iters = []
        setup_iter = []

        loop_end = []
        iter_filter = []

        for v in range(0, len(valid_dims)):
            dim, iters = valid_dims[v]
            limits = iters[0].limits
            symbols = flatten([i.expr_symbols for i in iters])

            has_sub_block = v < len(valid_dims) - 1

            l_idx = "((threadIdx.x %s) %% _block_%s)" % ("" if v == len(valid_dims) - 1 else ("/ (%s)" % ' * '.join("_block_%s" % x for x in dim_vars[v+1:len(valid_dims)])), dim_vars[v])
            kernel.append(c.Initializer(c.Value('int', dim.name + ("_0" if has_sub_block else "")), "blockIdx.%s * _block_%s %s+ %s" % (dim_vars[v], dim_vars[v], "* _sub_block_" + dim_vars[v] + " " if has_sub_block else "", l_idx)))          
            args = args.union(symbols)

            if has_sub_block:
                sub_var = "_sub_block_%s" % dim_vars[v]
                sub_iterator = "_" + dim_vars[v] + dim_vars[v]
                setup_iter.append(c.Initializer(c.Value("int", dim.name), "%s + %s" % (dim.name + "_0", sub_iterator)))

            # Add the iteration conditions
            iter_filter.append(c.If("%s < %s || %s > %s" % (dim.name, str(limits[0]), dim.name, str(limits[1])), c.Statement("continue") if v < len(valid_dims) - 1 else c.Statement("return")))

        body = setup_iter + iter_filter + [iet]

        for v in reversed(range(0, len(sub_blocks))):
            dim, iters = valid_dims[v]

            sub_var = "_sub_block_%s" % dim_vars[v]
            sub_iterator = "_" + dim_vars[v] + dim_vars[v]
            body = [c.Line("#pragma unroll"), c.Line("for (int %s = 0; %s < %s; %s++) {" % (sub_iterator, sub_iterator, sub_var, sub_iterator))] + body + [c.Line("}")]

        # todo: figure out something better based on looking at access for spatial reuse
        block = [1] * len(valid_dims)
        block[-1] = 32 # always want at least one warp worth, and preferably a multiple of warps
        if len(block) == 3:
            block[0] = 4
            block[1] = 4
        elif len(block) == 2:
            block[0] = 16
        else:
            block[0] = 128


        # Add the iteration body
        kernel.extend(as_tuple(body))

        # 'preferred' block is so named because we have no idea until at runtime how many
        # registers the CUDA compiler will use and thus the range of valid block sizes
        cuda_callable = CudaCallable(name=name, body=kernel, parameters=args, defines=[x[0] for x in valid_dims], preferred_block=block, preferred_sub_block=sub_blocks)

        return (cuda_callable, list([x[1][0] for x in valid_dims]))

    def _make_reductions(self, partree, force_atomic=False):
        if not force_atomic and not any(i.is_ParallelAtomic for i in FindNodes(Iteration).visit(partree)):
            return partree

        exprs = [i for i in FindNodes(Expression).visit(partree) if i.is_reduction]
        reductions = [(i.output, i.operation) for i in exprs]

        test0 = all(not i.is_Indexed for i, _ in reductions)

        if test0:
            # Implement reduction
            mapper = {partree.root: partree.root._rebuild(reduction=reductions)}
        elif all(i is OpInc for _, i in reductions):
            # Use atomic increments
            mapper = {i: i._rebuild(atomic=True) for i in exprs}
        else:
            raise NotImplementedError

        partree = Transformer(mapper).visit(partree)

        return partree

def tuple_to_dim3(grid):
    return "dim3(%s)" % ', '.join(str(x) if '/' not in str(x) else ("max(1, %s)" % str(x)) for x in grid )

@iet_pass
def kernel_tuning(iet):
    if not isinstance(iet, EntryFunction):
        return iet, {}
    kernel_calls = FindNodes(CudaCall).visit(iet.body)

    grouped = {n : list(set([c for c in kernel_calls if c.name == n])) for n in set([x.name for x in kernel_calls]) }
    unique = []
    non_unique = []
    for k, v in grouped.items():
        if len(v) == 1 or len(set([x.template_arguments for x in v])) == 1:
            unique += v
        else:
            non_unique += v
    tunes = []
    unique = sorted(unique, key=lambda x: x.name)
    non_unique = sorted(non_unique, key=lambda x: x.name)

    for call in unique:
        if call.preferred_block is None:
            continue
        setup_lambda = "[&](dim3 block, dim3 sub_block) {\n"
        suffix = ['.x', '.y', '.z']
        block_parameters = ['block' + suffix[i] for i in range(len(call.preferred_block))]
        subblock_parameters = ['sub_block' + suffix[i] for i in range(len(call.preferred_sub_block) if call.preferred_sub_block is not None else 0)]
        setup_lambda += 'return program.kernel("%s").instantiate(%s); }' % (call.name, ','.join(block_parameters + subblock_parameters + [ccode(x.rhs) for x in call.template_arguments]))
        preferred_sub_block = call.preferred_sub_block or []
        tunes.append(c.Line("""auto %s_tune = performTuning(_kernelTuning, "%s", %s, %s, %s, %d, %s);""" % (call.name, call.name, 
                                                                                                        tuple_to_dim3(call.preferred_block),
                                                                                                        tuple_to_dim3(preferred_sub_block),
                                                                                                        tuple_to_dim3(call.grid),
                                                                                                        len(call.preferred_block),
                                                                                                        setup_lambda))
                    )
    
    for call in non_unique:
        preferred_sub_block = call.preferred_sub_block or []
        tunes.append(c.Line("""auto %s_tune = std::make_pair(dim3(%s), dim3(%s))));""" % (call.name, 
                                                                                          ','.join([str(x) for x in call.preferred_block]),
                                                                                          ','.join([str(x) for x in preferred_sub_block]))))
    return iet._rebuild(body=iet.body._rebuild(body=flatten(tunes+[iet.body.body]))), {}

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
        self.gpu_nofit = options['gpu-nofit']

    def _alloc_local_array_on_high_bw_mem(self, site, obj, storage, devicerm=None):
        """
        Allocate a local Array in the device high bandwidth memory.
        """
        # Create a local static pointer so we can persist between runs
        decl = Definition(obj, initvalue="nullptr")
        doalloc = self.lang['device-alloc']
        dofree = self.lang['device-free']

        nbytes = SizeOf(obj._C_typedata)*obj.size
        init = doalloc(nbytes, None, retobj=obj._C_symbol)
        #allocs = (decl, Conditional(CondEq(obj._C_symbol, NullPointer()), init))
        allocs = (Call("PER_DEVICE_TEMP_GET", (ReservedWord(str(obj._C_typedata)), obj._C_symbol, nbytes)))

        #free = dofree(obj._C_name, None)
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

        nbytes_arg = SizeOf(obj.indexed._C_typedata)*obj.size

        alloc = List(body=[
            Call("PER_DEVICE_ARRAY_TEMP_DECLARE", (obj._C_symbol, ReservedWord(obj._C_typedata))),
            Call("PER_DEVICE_ARRAY_TEMP_GET", (obj._C_symbol, nbytes_arg), retobj=obj),
        ])

        free = Conditional(DeviceRM(), Call("PER_DEVICE_ARRAY_TEMP_DESTROY", (obj._C_symbol,)))

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

    @iet_visit
    def derive_cuda_casts(self, iet, **kwargs):
        # Don't generate unnecessary casts in CUDA kernels
        kernels = FindNodes(CudaCallable).visit(iet)
        calls = FindNodes(CudaCall).visit(iet)
        mapper = {}
        
        for kernel in kernels:
            indexeds = FindSymbols('indexeds|indexedbases').visit(kernel)
            defines = set(FindSymbols('defines').visit(kernel)) - set(kernel.parameters)
            bases = sorted({i.base for i in indexeds}, key=lambda i: i.name)
            casts = [CudaKernelPointerCast(i.function, obj=i) for i in bases
                    if i.function not in defines]

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
            cuda_filter = lambda n: isinstance(n, CudaCall) or isinstance(n, CudaCallable) or isinstance(n, CudaDealloc) or isinstance(n, PragmaTransfer) or isinstance(n, CudaHostFuncCall) or isinstance(n, CudaTransfer)
            # Candidates
            indexeds = FindSymbols('indexeds|indexedbases', stop_filter=cuda_filter).visit(iet)

            # Create Function -> n-dimensional array casts
            # E.g. `float (*u)[.] = (float (*)[.]) u_vec->data`
            # NOTE: a cast is needed only if the underlying data object isn't already
            # defined inside the kernel, which happens, for example, when:
            # (i) Dereferencing a PointerArray, e.g., `float (*r0)[.] = (float(*)[.]) pr0[.]`
            # (ii) Declaring a raw pointer, e.g., `float * r0 = NULL; *malloc(&(r0), ...)
            # we use iet.body here because we manually futz with the defines for some nodes to
            # coerce Devito into outputting function signatures the way we want them
            defines = set(FindSymbols('defines', stop_filter=cuda_filter).visit(iet.body))
            bases = sorted({i.base for i in indexeds}, key=lambda i: i.name)
            casts = [self.lang.PointerCast(i.function, obj=i) for i in bases
                    if i.function not in defines]

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

    @iet_pass
    def place_cuda_casts(self, iet, **kwargs):
        return Transformer(kwargs['mapper']).visit(iet), {}
    
    def process(self, graph):
        """
        Apply the `place_transfers`, `place_definitions` and `place_casts` passes.
        """
        mapper = self.derive_transfers(graph)
        self.place_transfers(graph, mapper=mapper)
        self.place_definitions(graph)
        cast_mapper = self.derive_cuda_casts(graph)
        cast_mapper = flatten_dict(cast_mapper, prefix=None)
        cast_mapper = flatten_dict(self.derive_cuda_kernel_call_parameters(graph, mapper=cast_mapper), prefix=None)
        self.place_cuda_casts(graph, mapper=cast_mapper)
        self.place_cuda_non_kernel_casts(graph)
        kernel_tuning(graph)
        
        self.tidy_up(graph)

def realign_iet(iet):
    """
    Attempt to realign the final iteration dimensions in an IET to better suit the GPU's cachelines.

    eg. 
    for (xi = x_m; xi < x_M; xi++)
    for (yi = y_m; yi < y_M; yi++)
        d[xi + 8][yi + 7] = sqrt(d[xi + 8][yi + 6])
        d2[xi + 8][yi + 7] = sqrt(d2[xi + 8][yi + 7])

    becomes

    for (xi = x_m + 8; xi < x_M + 8; xi++)
    for (yi = y_m + 7; yi < y_M + 7; yi++)
        d[xi + 8][yi] = sqrt(d[xi + 8][yi - 1])
        d2[xi + 8][yi] = sqrt(d2[xi + 8][yi])

    This way, the CUDA threads in each warp are reading and writing a centre point that is
    aligned with the 256-byte aligned grid. Granted, they'll usually also be reading misaligned
    points along the most frequently-changing axis in realistic workloads, 
    """

    # Find the functions we're writing to
    exs = FindNodes(Expression).visit(iet)
    exprs = [x for x in exs if x.output.is_Indexed]

    # Find all the dimensions they write to
    all_dims = set(flatten([expr.expr.lhs.function.dimensions[-1] for expr in exs if len(expr.expr.lhs.function.dimensions) > 0]))
    index_dims = set(flatten([[d for d in expr.expr.dimensions] for expr in exs]))
    index_map = {d.root : d for d in index_dims}
    # and reference those to the indexed dimensions
    used_index_map = {d: list(set(flatten([expr.expr.lhs.indices[d] for expr in exprs if d in expr.expr.lhs.indices._getters]))) for d in index_map.keys()}

    # for now, give up early if anything has multiple used indices
    if any([len(x) > 1 for x in used_index_map.values()]):
        return iet

    offset_map = {d: (-uxreplace(used_index_map[d][0], {index_map[d]: Number(0) })) for d in all_dims if d in used_index_map and len(used_index_map[d]) > 0}
    ioffset_map = {index_map[d]: offset_map[d] for d in offset_map }

    # Rebuild the expressions, and include the numeric offsets into the iterations
    # The simplify() calls are just to make the generated code a little more readable.
    # May want to remove them if they're taking an undue amount of time to process 
    # (or find a way of having it just simplify the surrounding terms?)

    return IterationLimitTranslator(ioffset_map).visit(iet)

class IterationLimitTranslator(Visitor):
    def __init__(self, dim_mapper):
        super(Visitor, self).__init__()
        self._dim_mapper = dim_mapper
        self._expr_map = { k: k + v for k, v in dim_mapper.items() }

    def visit_object(self, o, **kwargs):
        return o

    def visit_tuple(self, o, **kwargs):
        visited = tuple(self._visit(i, **kwargs) for i in o)
        return tuple(i for i in visited if i is not None)

    visit_list = visit_tuple

    def visit_Iteration(self, o, **kwargs):
        if o.dim in self._dim_mapper:
            return o._rebuild(limits=(o.limits[0] - self._dim_mapper[o.dim],
                                      o.limits[1] - self._dim_mapper[o.dim],
                                      1),
                              nodes=self._visit(o.nodes, **kwargs))

        else:
            children = [self._visit(i, **kwargs) for i in o.children]
            return o._rebuild(*children, **o.args_frozen)

    def visit_Expression(self, o, **kwargs):
        return o._rebuild(expr=o.expr.func(
            simplify(uxreplace(o.expr.lhs, self._expr_map)), 
            pow_to_mul(simplify(uxreplace(o.expr.rhs, self._expr_map))),
            ispace=o.expr.ispace.translate(self._dim_mapper))
        )
        
    def visit_Node(self, o, **kwargs):
        children = [self._visit(i, **kwargs) for i in o.children]
        return o._rebuild(*children, **o.args_frozen)

def flatten_dict(dd, separator='_', prefix=''):
    return { k : v
             for kk, vv in dd.items()
             for k, v in flatten_dict(vv, separator, kk).items()
             } if isinstance(dd, dict) else { prefix : dd }

@iet_pass
def preload_kernels(iet):
    calls = FindNodes(CudaCall).visit(iet)
    unique_kernels = set([(c.name, c.template_arguments) for c in calls])

class CudaHostFuncCall(AsyncCall):
    def __init__(self, name, arguments=None, retobj=None, is_indirect=False,
                 cast=False, writes=None, types=None, declares=True, stream=None):
        super().__init__(name, arguments=arguments, retobj=retobj, is_indirect=is_indirect,
                         cast=cast, writes=writes, types=types, declares=declares)
        self.stream = stream

class CudaHostFuncCallable(AsyncCallable):
    pass

class CudaOrchestrator(Orchestrator):
    lang = CudaBB

    _memcpy_stream = MemCopyStream()
    _host_stream = HostStream()
    _kernel_stream = KernelStream()

    def _make_waitevent(self, iet, sync_ops):
        waitloop = List(
            header=c.Comment("Wait for `%s` to be copied to the host" %
                             ",".join(s.function.name for s in sync_ops)),
            body=[self.lang._map_wait_recreate_event(s, stream=self._kernel_stream) for s in sync_ops],
            footer=c.Line()
        )

        iet = List(body=(waitloop,) + iet.body)

        return iet, []

    def _make_withlock(self, iet, sync_ops):
        preactions = [c.Comment("Block the copy until it's safe"),
                      BlankLine]

        # the main kernel stream should mark this as the appropriate place for it to start
        preactions.extend([self.lang._map_fire_event(s, stream=self._kernel_stream) for s in sync_ops])
        # these should run on the memcpy stream, so it needs to wait
        preactions.extend([self.lang._map_wait_event(s, stream=self._memcpy_stream) for s in sync_ops])

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
                iet = iet._rebuild(body=List(body=[iet.body, CudaChecked(Call("cudaDeviceSynchronize", None))]))
            return iet, {}

        callbacks = OrderedDict([
            (WithLock, self._make_withlock),
            (WaitLock, self._make_waitlock),
            (FetchUpdate, self._make_fetchupdate),
            (PrefetchUpdate, self._make_prefetchupdate),
            (ReleaseLock, self._nop)
        ])

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

        events = [List(body=[Definition(e, None, None, NullPointer()),
                             self.lang.mapper['create-event'](e._C_symbol),
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

        initialization.append(Definition(sdata, None, None, NullPointer()))

        # Activation
        #if threads.size == 1:
        d = 0

        activation = [c.Comment("Allocate a new block of data for this invocation"),
                      Call("posix_memalign", (VOID(Byref(sdata), '**'), 64, SizeOf(sdata._C_typedata))),
                      call0]
        activation.extend([DummyExpr(FieldFromComposite(i.name, sdata[d]), i)
                           for i in sdata.ncfields])

        activation.append(
           c.Statement("cudaLaunchHostFunc(%s, (cudaHostFn_t)%s, %s)" % (n.stream if n.stream is not None else 0, n.name, ccode(sbase + d))),
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
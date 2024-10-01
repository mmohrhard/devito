from collections import OrderedDict
from devito.cuda.visitors import IterationExtractor
from devito.ir.equations.equation import OpInc
from devito.ir.iet.efunc import EntryFunction
from devito.ir.iet.nodes import Callable, Expression, Iteration
from devito.ir.iet.utils import filter_iterations, retrieve_iteration_tree
from devito.ir.iet.visitors import FindNodes, FindSymbols, Transformer
from devito.logger import debug
from devito.passes.iet.misc import is_on_device
from devito.passes.iet.parpragma import PragmaShmTransformer


from devito.sycl.lang import SyclBB
from devito.sycl.nodes import SyclAtomicExpression, SyclKernel, SyclKernelLaunch
from devito.tools.utils import as_tuple, flatten


import cgen as c

__all__ = ["DeviceSyclizer"]


class DeviceSyclizer(PragmaShmTransformer):
    lang = SyclBB
    DeviceIteration = lang.DeviceIteration

    count = 0

    def __init__(self, sregistry, options, platform, compiler):
        super().__init__(sregistry, options, platform, compiler)

        self.gpu_fit = options["gpu-fit"]
        self.gpu_nofit = options["gpu-nofit"]

    def _extract_kernels(self, candidates, nthreads=None):
        assert candidates

        root = candidates[0]
        if self._is_offloadable(root):
            kernel_name = f"{self.kernel_basename}{self.count}"

            kernel, extracted_iterators = self._make_sycl_kernel(kernel_name, root)

            # find the non-derived dimensions we're iterating over, since the dimension
            # list for an Iteration includes the original dimension and the
            # derived version
            kdims = [
                next(filter(lambda x: x.is_Derived is False, c.dimensions)).symbolic_size
                for c in extracted_iterators
            ][:3]

            # If we don't find any non-derived dimensions, use the derived ones I guess?
            if len(kdims) == 0:
                kdims = [c.dimensions[0].symbolic_size for c in extracted_iterators][:3]

            if len(kdims) == 0:
                return root, None, None

            # The GPU wants to iterate over the whole problem space for
            # memory alignment reasons;
            # filtering out unwanted points in the space on the GPU
            # is very approximately zero-cost

            partree = SyclKernelLaunch(name=kernel_name, kernel=kernel, dims=kdims)

            # partree = CudaCall(
            #     kernel_name,
            #     kgrid,
            #     kthread,
            #     preferred_block=kernel.preferred_block,
            #     preferred_sub_block=kernel.preferred_sub_block,
            #     arguments=kernel.parameters,
            #     kernel=kernel,
            #     stream=KernelStream(),
            # )
            # Make sure that the enclosing function knows we need the full size
            # of the Functions
            # partree.expr_symbols = as_tuple(flatten((partree.expr_symbols, kdims)))

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

        debug("in syclizer")
        # Name kernels according to the name of the EntryFunction by default so that
        # profiling multiple operators in a single Nsight run produces more
        # meaningful summary data
        try:
            self.kernel_basename = FindNodes(EntryFunction).visit(iet)[0].name + "_kernel"
        except IndexError:
            self.kernel_basename = "kernel"

        for tree in retrieve_iteration_tree(iet, mode="superset"):
            # Get the parallelizable Iterations in `tree`
            candidates = filter_iterations(tree, key=self.key)
            if not candidates:
                continue

            # Outer parallelism
            root, partree, _ = self._extract_kernels(candidates)
            if partree is None or root in mapper:
                continue

            mapper[root] = partree

        iet = Transformer(mapper).visit(iet)

        attrs = {
            "efuncs": kernels,
            "includes": self.lang["headers"],
            "globals": self.lang["global-decls"],
        }

        return iet, attrs

    def _make_nested_partree(self, partree):
        if isinstance(partree, Callable) or isinstance(
            partree.root, self.DeviceIteration
        ):
            # no-op for now
            return partree
        else:
            return super()._make_nested_partree(partree)

    def _make_reductions(self, partree, force_atomic=False):
        if not force_atomic and not any(
            i.is_ParallelAtomic for i in FindNodes(Iteration).visit(partree)
        ):
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

    def _make_sycl_kernel(self, name, body):
        # Find the iterators we consider eligible for being the GPU grid dimensions
        iterations = list(
            [i for i in FindNodes(Iteration).visit(body) if i.is_ParallelRelaxed]
        )
        # FIXME: string matching is bad; use properties?
        possible_iter_dimensions = list(
            OrderedDict.fromkeys(
                [x.dim for x in iterations if not x.dim.name.startswith("par_dim")]
            )
        )

        grouped_iters = [
            (x, list(OrderedDict.fromkeys([i for i in iterations if i.dim == x])))
            for x in possible_iter_dimensions
        ]

        valid_dims = list(
            filter(lambda i: len(set([z.limits for z in i[1]])) == 1, grouped_iters)
        )
        if len(valid_dims) > 3:
            valid_dims = valid_dims[0:3]

        iet = IterationExtractor([d[0] for d in valid_dims]).visit(body)

        # if any of the iterations we're extracting have the atomic flag set,
        # we need to force all reductions in the kernel to be atomic regardless
        # of any inner iterations
        force_atomic = any(
            i.is_ParallelAtomic for i in flatten([d[1] for d in valid_dims])
        )
        iet = self._make_reductions(iet, force_atomic=force_atomic)
        # replace any atomic ops
        exprs = [e for e in FindNodes(Expression).visit(iet) if e.is_atomic]
        mapper = dict(
            [
                (
                    i,
                    SyclAtomicExpression(i.expr, i.pragmas, i.init, i.operation),
                )
                for i in exprs
            ]
        )
        iet = Transformer(mapper).visit(iet)

        kernel = []
        args = set()

        # TODO: tune sub-blocking
        sub_blocks = [1] * (len(valid_dims) - 1)
        if len(sub_blocks) > 0:
            sub_blocks[0] = 2

        # TODO: figure out something better based on looking at access for spatial reuse
        block = [1] * len(valid_dims)
        # always want at least one warp worth, and preferably a
        # multiple of warps
        block[-1] = 32

        if len(block) == 3:
            block[0] = 1
            block[1] = 32
        elif len(block) == 2:
            block[0] = 16
        else:
            block[0] = 128

        setup_iter = []

        iter_filter = []

        for v in range(0, len(valid_dims)):
            dim, iters = valid_dims[v]
            limits = iters[0].limits
            symbols = flatten([i.expr_symbols for i in iters])

            kernel.append(c.Initializer(c.Value("int", dim.name), "_item[%d]" % v))
            args = args.union(symbols)

            # Add the iteration conditions
            iter_filter.append(
                c.If(
                    "%s < %s || %s > %s"
                    % (dim.name, str(limits[0]), dim.name, str(limits[1])),
                    c.Statement("return"),
                )
            )

        body = setup_iter + iter_filter + [iet]

        # Add the iteration body
        kernel.extend(as_tuple(body))

        # 'preferred' block is so named because we have no idea until at runtime how many
        # registers the CUDA compiler will use and thus the range of valid block sizes
        sycl_kernel = SyclKernel(
            name=name, body=kernel, parameters=args, defines=[x[0] for x in valid_dims]
        )

        return (sycl_kernel, list([x[1][0] for x in valid_dims]))

    def _is_offloadable(self, iet):
        """
        True if the IET computation is offloadable to device, False otherwise.
        """
        expressions = FindNodes(Expression).visit(iet)
        if any(
            not is_on_device(e.write, self.gpu_fit, self.gpu_nofit) for e in expressions
        ):
            return False

        functions = FindSymbols().visit(iet)
        buffers = [f for f in functions if f.is_Array and f._mem_mapped]
        hostfuncs = [
            f for f in functions if not is_on_device(f, self.gpu_fit, self.gpu_nofit)
        ]
        return not (buffers and hostfuncs)

from collections import OrderedDict
from devito.passes.iet.parpragma import PragmaDeviceAwareTransformer
from devito.ir.iet import (
    FindNodes,
    EntryFunction,
    Transformer,
    Iteration,
    Expression,
    Callable,
)
from devito.ir.iet.utils import filter_iterations, retrieve_iteration_tree
from devito.ir.equations import OpInc
from devito.tools import as_tuple, flatten

from devito.cuda.nodes import CudaCall, KernelStream, CudaAtomicExpression, CudaCallable

from devito.cuda.passes.realign import realign_iet
from devito.cuda.lang import CudaBB
from devito.cuda.visitors import IterationExtractor

import cgen as c

__all__ = ["DeviceCudaizer"]


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
            kgrid = kdims.copy()

            # the grid/threads are (for now) set up in some C++ code from a header
            kthread = [1] * len(kdims)

            partree = CudaCall(
                kernel_name,
                kgrid,
                kthread,
                preferred_block=kernel.preferred_block,
                preferred_sub_block=kernel.preferred_sub_block,
                arguments=kernel.parameters,
                kernel=kernel,
                stream=KernelStream(),
            )
            # Make sure that the enclosing function knows we need the full size
            # of the Functions
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

        for tree in retrieve_iteration_tree(iet, mode="superset"):
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

    def _make_cuda_kernel(self, name, body):
        body = realign_iet(body)
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
                    CudaAtomicExpression(i.expr, i.pragmas, i.init, i.operation),
                )
                for i in exprs
            ]
        )
        iet = Transformer(mapper).visit(iet)

        # Now, generate the iteration dimension variables from the blockIdx/threadIdx
        # These end up reversed because warp thread order in CUDA for >1D is column-major
        # and we want contiguous memory access
        dim_vars = ["x", "y", "z"]

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

            has_sub_block = v < len(valid_dims) - 1

            # TODO: Don't just jam C++ in here; turn it into nodes that we lower
            # into C++ at codegen time
            l_idx = "((threadIdx.x %s) %% _block_%s)%s" % (
                (
                    ""
                    if v == len(valid_dims) - 1
                    else (
                        "/ (%s)"
                        % " * ".join(
                            "_block_%s" % x for x in dim_vars[v + 1 : len(valid_dims)]
                        )
                    )
                ),
                dim_vars[v],
                "* _sub_block_" + dim_vars[v] + " " if has_sub_block else "",
            )
            kernel.append(
                c.Initializer(
                    c.Value("int", dim.name + ("_0" if has_sub_block else "")),
                    "blockIdx.%s * _block_%s %s+ %s"
                    % (
                        dim_vars[v],
                        dim_vars[v],
                        ("* _sub_block_" + dim_vars[v] + " " if has_sub_block else ""),
                        l_idx,
                    ),
                )
            )
            args = args.union(symbols)

            if has_sub_block:
                sub_iterator = "_" + dim_vars[v] + dim_vars[v]
                setup_iter.append(
                    c.Initializer(
                        c.Value("int", dim.name),
                        "%s + %s" % (dim.name + "_0", sub_iterator),
                    )
                )

            # Add the iteration conditions
            iter_filter.append(
                c.If(
                    "%s < %s || %s > %s"
                    % (dim.name, str(limits[0]), dim.name, str(limits[1])),
                    c.Statement("continue") if has_sub_block else c.Statement("return"),
                )
            )

        body = setup_iter + iter_filter + [iet]

        for v in reversed(range(0, len(sub_blocks))):
            dim, iters = valid_dims[v]

            sub_var = "_sub_block_%s" % dim_vars[v]
            sub_iterator = "_" + dim_vars[v] + dim_vars[v]
            body = (
                [
                    c.Line("#pragma unroll"),
                    c.Line(
                        "for (int %s = 0; %s < %s; %s++) {"
                        % (sub_iterator, sub_iterator, sub_var, sub_iterator)
                    ),
                ]
                + body
                + [c.Line("}")]
            )

        # Add the iteration body
        kernel.extend(as_tuple(body))

        # 'preferred' block is so named because we have no idea until at runtime how many
        # registers the CUDA compiler will use and thus the range of valid block sizes
        cuda_callable = CudaCallable(
            name=name,
            body=kernel,
            parameters=args,
            defines=[x[0] for x in valid_dims],
            preferred_block=block,
            preferred_sub_block=sub_blocks,
        )

        return (cuda_callable, list([x[1][0] for x in valid_dims]))

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

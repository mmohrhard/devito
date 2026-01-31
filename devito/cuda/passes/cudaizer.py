from collections import OrderedDict

import cgen as c

from devito.cuda.lang import CudaBB
from devito.cuda.nodes import (
    CudaAtomicExpression,
    CudaCall,
    CudaCallable,
    KernelStream,
)
from devito.cuda.passes.realign import realign_iet
from devito.cuda.visitors import IterationExtractor
from devito.ir.equations import OpInc
from devito.ir.iet import (
    Callable,
    EntryFunction,
    Expression,
    FindNodes,
    Iteration,
    Transformer,
)
from devito.ir.iet.utils import filter_iterations, retrieve_iteration_tree
from devito.passes.iet.parpragma import PragmaDeviceAwareTransformer
from devito.tools import as_tuple, flatten

__all__ = ["DeviceCudaizer"]


class DeviceCudaizer(PragmaDeviceAwareTransformer):
    lang = CudaBB
    DeviceIteration = lang.DeviceIteration

    count = 0

    def __init__(self, sregistry, options, platform, compiler):
        super().__init__(sregistry, options, platform, compiler)
        self._realign_iet_opt = options.get("realign_iet", True)
        self._unroll_sub_blocks = options.get("cuda-unroll-subblocks", True)
        self._block_sizes = options.get("cuda-par-block-sizes", [])

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
                next(filter(lambda x: x.is_Derived is False, c.dimensions))
                for c in extracted_iterators
            ][:3]

            # If we don't find any non-derived dimensions, use the derived ones I guess?
            if len(kdims) == 0:
                kdims = [c.dimensions[0] for c in extracted_iterators][:3]

            if len(kdims) == 0:
                return root, None, None

            # The GPU wants to iterate over the whole problem space for
            # memory alignment reasons;
            # filtering out unwanted points in the space on the GPU
            # is very approximately zero-cost
            #
            # justinw 12/25: Oops, we found an edge case! MPI halo exchanges
            # call into the overlapping kernel 26 times with small iteration spaces
            # (6 faces + 8 corners + 12 edges), and for any reasonably sized grid
            # we end up with substantial kernel launch overhead.
            #
            # We now pass the dimension symbolic minimums/maximums to the launch
            # function, and it'll round those down to the nearest block size, and
            # launch the smallest possible grid
            kgrid = [d.symbolic_size for d in kdims].copy()
            kmins = [d.symbolic_min for d in kdims].copy()
            kmaxs = [d.symbolic_max for d in kdims].copy()

            # the grid/threads are (for now) set up in some C++ code from a header
            kthread = [1] * len(kdims)

            partree = CudaCall(
                kernel_name,
                kgrid,
                kmins,
                kmaxs,
                kthread,
                preferred_block=kernel.preferred_block,
                preferred_sub_block=kernel.preferred_sub_block,
                arguments=kernel.parameters,
                kernel=kernel,
                stream=KernelStream(),
            )

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
        if self._realign_iet_opt is True:
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
            (
                x,
                list(OrderedDict.fromkeys([i for i in iterations if i.dim == x])),
            )
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
        dim_vars = ["z", "y", "x"]

        kernel = []
        args = set()

        # TODO: tune sub-blocking
        sub_blocks = [1, 1, 1]
        num_subblocks = max(0, len(valid_dims) - 1)

        # TODO: figure out something better based on looking at access for spatial reuse
        block = [1, 1, 1]
        # always want at least one warp worth, and preferably a
        # multiple of warps
        block[0] = 32

        if len(valid_dims) == 3:
            block[2] = 2
            block[1] = 4
            sub_blocks[0] = 1
        elif len(valid_dims) == 1:
            block[0] = 128

        if not any(s > 1 for s in sub_blocks):
            num_subblocks = 0

        # allow overriding in options
        if isinstance(self._block_sizes, dict):
            # filtered by kernel name
            for n in self._block_sizes:
                if n in name:
                    bs = self._block_sizes[n]
                    if (
                        len(bs) == 2
                        and isinstance(bs[0], tuple)
                        and isinstance(bs[1], tuple)
                    ):
                        for d in range(0, 3):
                            if len(bs[0]) > d:
                                block[d] = bs[0][d]
                            else:
                                block[d] = 1

                            if len(bs[1]) > d:
                                sub_blocks[d] = bs[1][d]
                            else:
                                sub_blocks[d] = 1

                            if sub_blocks[d] > 1:
                                num_subblocks = max(num_subblocks, 3 - d)
        elif (
            isinstance(self._block_sizes, tuple)
            and len(self._block_sizes) == 2
            and isinstance(self._block_sizes[0], tuple)
            and isinstance(self._block_sizes[1], tuple)
        ):
            # generic for all kernels
            bs = self._block_sizes
            for d in range(0, 3):
                if len(bs[0]) > d:
                    block[d] = bs[0][d]
                else:
                    block[d] = 1

                if len(bs[1]) > d:
                    sub_blocks[d] = bs[1][d]
                else:
                    sub_blocks[d] = 1

                if sub_blocks[d] > 1:
                    num_subblocks = max(num_subblocks, 3 - d)

        for idx in range(0, 3):
            if sub_blocks[-(1 + idx)] > 1:
                num_subblocks = 3 - idx
                break

        num_subblocks = min(len(valid_dims), num_subblocks)
        setup_iter = []

        iter_filter = []

        for v in range(0, len(valid_dims)):
            dv = v - len(valid_dims)
            dim, iters = valid_dims[v]
            limits = iters[0].limits
            symbols = flatten([i.expr_symbols for i in iters])

            has_sub_block = v < num_subblocks

            # TODO: Don't just jam C++ in here; turn it into nodes that we lower
            # into C++ at codegen time
            l_idx = "((threadIdx.x%s) %% _block_%s)%s" % (
                (
                    ""
                    if v == len(valid_dims) - 1
                    else (
                        " / (%s)"
                        % " * ".join("_block_%s" % x for x in dim_vars[dv + 1 :])
                    )
                ),
                dim_vars[dv],
                " * _sub_block_" + dim_vars[dv] + " " if has_sub_block else "",
            )
            kernel.append(
                c.Initializer(
                    c.Value("int", dim.name + ("_0" if has_sub_block else "")),
                    "(blockIdx.%s * _block_%s %s) + %s + %s"
                    % (
                        dim_vars[dv],
                        dim_vars[dv],
                        (" * _sub_block_" + dim_vars[dv] + " " if has_sub_block else ""),
                        "_offsets." + dim_vars[v],
                        l_idx,
                    ),
                )
            )
            args = args.union(symbols)

            if has_sub_block:
                sub_iterator = "_" + dim_vars[dv] + dim_vars[dv]
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

        unroll_sub_blocks = self._unroll_sub_blocks
        if isinstance(self._unroll_sub_blocks, dict):
            for n in self._unroll_sub_blocks:
                if n in name:
                    unroll_sub_blocks = self._unroll_sub_blocks[n]

        for v in range(0, num_subblocks):
            dv = v - len(valid_dims)
            dim, iters = valid_dims[v]

            should_unroll = unroll_sub_blocks is True or (
                unroll_sub_blocks >= 1 and unroll_sub_blocks < (len(valid_dims) - v)
            )

            sub_var = "_sub_block_%s" % dim_vars[dv]
            sub_iterator = "_" + dim_vars[dv] + dim_vars[dv]
            body = (
                [
                    c.Line("#pragma unroll") if should_unroll else c.Line(""),
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
            preferred_block=block[0 : len(valid_dims)],
            preferred_sub_block=sub_blocks[0 : len(valid_dims)],
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

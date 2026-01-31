from typing import Any, Iterable

import cgen as c

from devito.cuda.nodes import CudaCall, CudaCallable
from devito.cuda.types import (
    JitifyKernelInstantiation,
    JitKernelTuningParams,
    JitOnly,
)
from devito.cuda.utils import tuple_to_dim3
from devito.ir.iet.efunc import EntryFunction
from devito.ir.iet.nodes import Callable, Definition, Node
from devito.ir.iet.visitors import FindNodes, Transformer
from devito.passes.iet.engine import Graph, iet_pass
from devito.symbolics.printer import ccode
from devito.tools.utils import as_tuple, flatten

__all__ = ["kernel_tuning"]


class CallTracker:
    def __init__(self):
        self.calls = []

    def add(self, o: CudaCall | Iterable[CudaCall]):
        if isinstance(o, Iterable):
            self.calls.extend(o)
        else:
            self.calls.append(o)


def kernel_tuning(graph: Graph, **kwargs):
    tracker = CallTracker()

    _gather_kernels(graph, tracker=tracker)

    all_params: dict[str, JitKernelTuningParams] = {}
    all_instantiations: dict[str, JitifyKernelInstantiation] = {}

    grouped = {
        n: list(set([c for c in tracker.calls if c.name == n]))
        for n in set([x.name for x in tracker.calls])
    }

    unique = []
    non_unique = []
    for k, v in grouped.items():
        if len(v) == 1 or len(set(flatten([x.template_arguments for x in v]))) == 1:
            unique += v
        else:
            non_unique += v

    # this has held true for all our current operators?
    if len(non_unique) > 0:
        print("found non-unique calls?")

    tunes: list[Any] = []
    unique = sorted(unique, key=lambda x: x.name)
    non_unique = sorted(non_unique, key=lambda x: x.name)

    for call in unique:
        if call.preferred_block is None:
            continue
        setup_lambda = "[&](dim3 block, dim3 sub_block) {\n"
        suffix = [".x", ".y", ".z"]
        block_parameters = ["block" + suffix[i] for i in range(len(call.preferred_block))]
        subblock_parameters = [
            "sub_block" + suffix[i]
            for i in range(
                len(call.preferred_sub_block)
                if call.preferred_sub_block is not None
                else 0
            )
        ]
        setup_lambda += '\treturn program.kernel("%s").instantiate(%s);\n\t}' % (
            call.name,
            ",".join(
                block_parameters
                + subblock_parameters
                + [
                    ("std::max<int>(1, %s)" % ccode(x.rhs))
                    for x in call.template_arguments
                ]
            ),
        )
        preferred_sub_block = call.preferred_sub_block or []

        builder_name = "_build_" + call.name
        tune_name = call.name + "_tune"
        tuned_name = call.name + "_tuned"

        params = JitKernelTuningParams(tune_name)
        instantiated = JitifyKernelInstantiation(tuned_name)
        all_params[call.name] = params
        all_instantiations[call.name] = instantiated
        tunes.extend(
            [
                c.Line(
                    "std::function<jitify::KernelInstantiation(dim3, dim3)> %s = %s;"
                    % (builder_name, setup_lambda)
                ),
                Definition(
                    params,
                    initvalue="devito::cuda::performTuning(_kernelTuning"
                    + ', "%s", %s, %s, %s, %d, %s);'
                    % (
                        call.name,
                        tuple_to_dim3(call.preferred_block),
                        tuple_to_dim3(preferred_sub_block),
                        tuple_to_dim3(call.grid),
                        len(call.preferred_block),
                        builder_name,
                    ),
                ),
                Definition(
                    instantiated,
                    initvalue="%s(%s, %s)"
                    % (
                        builder_name,
                        "std::get<0>(%s)" % tune_name,
                        "std::get<1>(%s)" % tune_name,
                    ),
                ),
                c.Line(),
            ]
        )

    tunes.append(c.If(JitOnly(), c.Statement("return 0")))

    _apply_kernel_tuning(
        graph,
        tunes=tunes,
        params=all_params,
        instantiations=all_instantiations,
        **kwargs,
    )


@iet_pass
def _gather_kernels(iet, tracker: CallTracker, **kwargs) -> tuple[Node, dict[str, Node]]:
    if not isinstance(iet, Callable):
        return iet, {}

    if isinstance(iet, CudaCallable):
        return iet, {}

    calls = FindNodes(CudaCall).visit(iet.body)
    tracker.add(calls)

    return iet, {}


@iet_pass
def _apply_kernel_tuning(
    iet,
    tunes: list[Any] | None = None,
    params: dict[str, JitKernelTuningParams] | None = None,
    instantiations: dict[str, JitifyKernelInstantiation] | None = None,
    **kwargs,
) -> tuple[Node, dict[str, Node]]:
    if not isinstance(iet, Callable):
        return iet, {}

    if isinstance(iet, CudaCallable):
        return iet, {}

    kernel_calls = FindNodes(CudaCall).visit(iet.body)

    tunes = tunes or []
    params = params or {}
    instantiations = instantiations or {}

    mapper = {}

    for original_call in kernel_calls:
        if original_call.name in params:
            mapper[original_call] = original_call._rebuild(
                jit_instantiation=instantiations[original_call.name],
                tune=params[original_call.name],
            )

    if isinstance(iet, EntryFunction):
        jitonly = [JitOnly()]
        iet = iet._rebuild(
            body=iet.body._rebuild(
                body=flatten(tunes + [Transformer(mapper).visit(iet.body.body)])
            ),
            parameters=as_tuple(flatten([iet.parameters] + jitonly)),
        )
    else:
        iet = iet._rebuild(
            body=iet.body._rebuild(body=flatten(Transformer(mapper).visit(iet.body.body)))
        )

    return iet, {}

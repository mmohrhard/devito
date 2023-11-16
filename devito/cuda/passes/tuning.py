import cgen as c

from devito.cuda.nodes import CudaCall
from devito.cuda.utils import tuple_to_dim3
from devito.ir.iet.efunc import EntryFunction
from devito.ir.iet.visitors import FindNodes
from devito.passes.iet.engine import iet_pass
from devito.symbolics.printer import ccode
from devito.tools.utils import flatten

__all__ = ["kernel_tuning"]


@iet_pass
def kernel_tuning(iet, **kwargs):
    if not isinstance(iet, EntryFunction):
        return iet, {}
    kernel_calls = FindNodes(CudaCall).visit(iet.body)

    grouped = {
        n: list(set([c for c in kernel_calls if c.name == n]))
        for n in set([x.name for x in kernel_calls])
    }

    unique = []
    non_unique = []
    for k, v in grouped.items():
        if len(v) == 1 or len(set(flatten([x.template_arguments for x in v]))) == 1:
            unique += v
        else:
            non_unique += v

    # this has held true for all our current operators?
    assert len(non_unique) == 0

    tunes = []
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

        tunes.extend(
            [
                c.Line(
                    "std::function<jitify::KernelInstantiation(dim3, dim3)> %s = %s;"
                    % (builder_name, setup_lambda)
                ),
                c.Line(
                    'auto %s = performTuning(_kernelTuning, "%s", %s, %s, %s, %d, %s);'
                    % (
                        tune_name,
                        call.name,
                        tuple_to_dim3(call.preferred_block),
                        tuple_to_dim3(preferred_sub_block),
                        tuple_to_dim3(call.grid),
                        len(call.preferred_block),
                        builder_name,
                    )
                ),
                c.Line(
                    "auto %s = %s(%s, %s);"
                    % (
                        tuned_name,
                        builder_name,
                        "std::get<0>(%s)" % tune_name,
                        "std::get<1>(%s)" % tune_name,
                    )
                ),
                c.Line(),
            ]
        )

    return (
        iet._rebuild(body=iet.body._rebuild(body=flatten(tunes + [iet.body.body]))),
        {},
    )

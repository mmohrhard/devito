import cgen as c

from devito.cuda.nodes import CudaCall
from devito.cuda.utils import tuple_to_dim3
from devito.ir.iet.efunc import EntryFunction
from devito.ir.iet.visitors import FindNodes
from devito.passes.iet.engine import iet_pass
from devito.symbolics.printer import ccode
from devito.tools.utils import flatten

__all__ = ['kernel_tuning']

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

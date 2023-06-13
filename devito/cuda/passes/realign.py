from sympy.core.numbers import Number
from sympy import simplify

from devito.ir.iet.nodes import Expression
from devito.ir.iet.visitors import FindNodes, Visitor
from devito.symbolics.manipulation import pow_to_mul, uxreplace
from devito.tools.utils import flatten

__all__ = ['realign_iet']

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
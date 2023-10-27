from devito.ir.iet.visitors import Visitor
from devito.ir.iet.nodes import List


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

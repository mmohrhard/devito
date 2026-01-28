from devito.ir.iet import TimedList
from devito.ir.iet.visitors import FindNodes, Transformer
from devito.cuda.nodes import CudaTimedList

__all__ = ["apply_cuda_timing"]


# Replace Devito TimedLists with CudaTimedLists that call into supporting
# code in the header to measure GPU kerneel execution time.
def apply_cuda_timing(iet):
    mapper = {}
    lists = FindNodes(TimedList).visit(iet)
    for tl in lists:
        mapper[tl] = CudaTimedList(tl.timer, tl.name, tl.body)

    return Transformer(mapper).visit(iet)

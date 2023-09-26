from devito.ir import FindSymbols, derive_parameters
from devito.tools import filter_ordered, flatten

__all__ = ["flatten_dict", "cuda_derive_parameters", "tuple_to_dim3"]


def flatten_dict(dd, separator="_", prefix=""):
    return (
        {
            k: v
            for kk, vv in dd.items()
            for k, v in flatten_dict(vv, separator, kk).items()
        }
        if isinstance(dd, dict)
        else {prefix: dd}
    )


def cuda_derive_parameters(iet):
    indexeds = FindSymbols("indexeds|indexedbases").visit(iet)
    bases = sorted({i.base for i in indexeds}, key=lambda i: i.name)

    return filter_ordered(
        flatten([derive_parameters(iet), FindSymbols("basics").visit(bases)])
    )


def tuple_to_dim3(grid):
    return "dim3(%s)" % ", ".join(
        str(x) if "/" not in str(x) else ("max(1, %s)" % str(x)) for x in grid
    )

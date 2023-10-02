from pathlib import Path

from devito.ir import FindSymbols, derive_parameters
from devito.tools import filter_ordered, flatten

from sympy import Mod

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


def get_header_include_path() -> Path:
    return Path(__file__).resolve().parent / "headers"


def q_has_modulo(expr):
    if isinstance(expr, Mod):
        return True
    elif expr.is_Atom:
        return False
    else:
        return any(q_has_modulo(a) for a in expr.args)

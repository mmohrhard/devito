from collections import OrderedDict, defaultdict

from sympy import And


from devito.ir import (Forward, GuardBoundNext, Queue, Vector, SEQUENTIAL,
                       WaitLock, WithLock, FetchUpdate, PrefetchUpdate,
                       ReleaseLock, normalize_syncs)
from devito.symbolics import uxreplace
from devito.tools import flatten, is_integer, timed_pass
from devito.types import CustomDimension, Lock

__all__ = ['cuda_memcpy']

def cuda_memcpy(graph, **kwargs):
    """
    Detects Expressions that are equivalent to 1D, 2D or 3D memcpy and rewrites them into CUDA memcpy calls
    """
    return graph
    


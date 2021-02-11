from collections import defaultdict

from sympy import Add, Mul, collect

from devito.passes.clusters.utils import cluster_pass
from devito.symbolics import estimate_cost, retrieve_scalars
from devito.tools import ReducerMap
from devito.ir.support import IterationSpace
from devito.ir.support.space import Interval
from devito.ir.clusters.cluster import Cluster
from devito.symbolics import uxreplace, xreplace_indices

__all__ = ['skewing']


"""
Minimum operation count of an expression so that aggressive factorization
is applied.
"""


@cluster_pass
def skewing(cluster, *args):

    """
    Skew the accesses along the time Dimension.
    Example:
    
    Transform
    
    for i = 2, n-1
        for j = 2, m-1
            a[i,j] = (a[a-1,j] + a[i,j-1] + a[i+1,j] + a[i,j+1]) / 4
        end for
    end for
    
    to
    
    for i = 2, n-1
        for j = 2+i, m-1+i
            a[i,j-i] = (a[a-1,j-i] + a[i,j-1-i] + a[i+1,j-i] + a[i,j+1-i]) / 4
        end for
    end for
    """
    processed = []
    
    
    itintervals = cluster.ispace.itintervals 
    itdimensions = cluster.ispace.itdimensions
    sub_iterators = cluster.ispace.sub_iterators
    directions = cluster.ispace.directions

        
    skew_dims = {i.dim for i in cluster.ispace.intervals if i.dim.is_Time}
    try:
        skew_dim = skew_dims.pop()
    except KeyError:
        # No time dimensions -> nothing to do
        return cluster
    if len(skew_dims) > 0:
        raise ValueError("More than 1 time dimensions. Aborting tt...")

    mapper, intervals = {}, []

    # Initializing a default time_dim index position in loop
    index = 0
    # Skew dim will not be none here:
    for i in cluster.ispace.itintervals:
        import pdb;pdb.set_trace()
        if i.dim.is_Time:
            intervals.append(Interval(i.dim, 0, 0))
            skew_dim = i.dim
            index = intervals.index(Interval(i.dim, 0, 0))
        elif index < cluster.ispace.itintervals.index(i):
            mapper[i.dim] = i.dim - skew_dim
            intervals.append(Interval(i.dim, skew_dim, skew_dim))
        else:
            intervals.append(i)

    processed = xreplace_indices(cluster.exprs, mapper)

    ispace = IterationSpace(intervals, sub_iterators, directions)
    cluster = Cluster(processed, ispace, cluster.dspace, guards=cluster.guards)

    return cluster.rebuild(processed)



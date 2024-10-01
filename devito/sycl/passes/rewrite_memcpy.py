__all__ = ["sycl_memcpy"]


def sycl_memcpy(graph, **kwargs):
    """
    Detects Expressions that are equivalent to 1D, 2D or 3D
    memcpy and rewrites them into SYCL memcpy calls
    """
    return graph

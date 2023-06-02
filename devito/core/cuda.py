from cached_property import cached_property
from devito.core.gpu import DeviceNoopOperator, DeviceAdvOperator, DeviceFsgOperator, DeviceCustomOperator
from devito.passes.iet import DeviceCudaTarget

__all__ = ['DeviceNoopCudaOperator','DeviceAdvCudaOperator', 'DeviceAdvCudaOperator', 'DeviceFsgCudaOperator',
           'DeviceCustomCudaOperator']
# CUDA

class DeviceCudaOperatorMixin(object):
    from devito.ir.iet.cuda_visitors import CudaCGen
    
    _Target = DeviceCudaTarget
    _CodeGen = CudaCGen
    @classmethod
    def _normalize_kwargs(cls, **kwargs):
        oo = kwargs['options']
        oo.pop('openmp', None)

        kwargs = super()._normalize_kwargs(**kwargs)
        oo['cuda'] = True

        return kwargs


class DeviceNoopCudaOperator(DeviceCudaOperatorMixin, DeviceNoopOperator):
    pass

class DeviceAdvCudaOperator(DeviceCudaOperatorMixin, DeviceAdvOperator):
    pass


class DeviceFsgCudaOperator(DeviceCudaOperatorMixin, DeviceFsgOperator):
    pass

class DeviceCustomCudaOperator(DeviceCudaOperatorMixin, DeviceCustomOperator):
    pass
    
    @classmethod
    def _make_iet_passes_mapper(cls, **kwargs):
        mapper = super()._make_iet_passes_mapper(**kwargs)
        mapper['cuda'] = mapper['parallel']
        mapper['pthreadify'] = mapper['cuda-events']
        mapper['linearize'] = mapper['cuda-linearize']
        return mapper

    @classmethod
    def _make_clusters_passes_mapper(cls, **kwargs):
        mapper = super()._make_clusters_passes_mapper(**kwargs)
        return mapper

    _known_passes = DeviceCustomOperator._known_passes + ('cuda', 'cuda-events', 'cuda-memcpy')
    assert not (set(_known_passes) & set(DeviceCustomOperator._known_passes_disabled))

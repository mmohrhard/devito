from devito.passes.iet import linearize, cuda_linearize
from devito.passes.iet.languages.C import CDataManager
from devito.passes.iet.languages.openmp import (SimdOmpizer, Ompizer, DeviceOmpizer,
                                                OmpDataManager, DeviceOmpDataManager,
                                                OmpOrchestrator)
from devito.passes.iet.languages.openacc import (DeviceAccizer, DeviceAccDataManager,
                                                 AccOrchestrator)
from devito.passes.iet.languages.cuda import (DeviceCudaizer, DeviceCudaDataManager, CudaOrchestrator)                                                 

__all__ = ['CTarget', 'OmpTarget', 'DeviceOmpTarget', 'DeviceAccTarget', 'DeviceCudaTarget']


class Target(object):
    Parizer = None
    DataManager = None
    Orchestrator = None
    Linearizer = linearize

class CTarget(Target):
    Parizer = SimdOmpizer
    DataManager = CDataManager
    Linearizer = linearize


class OmpTarget(Target):
    Parizer = Ompizer
    DataManager = OmpDataManager
    Linearizer = linearize

class DeviceOmpTarget(Target):
    Parizer = DeviceOmpizer
    DataManager = DeviceOmpDataManager
    Orchestrator = OmpOrchestrator
    Linearizer = linearize

class DeviceAccTarget(Target):
    Parizer = DeviceAccizer
    DataManager = DeviceAccDataManager
    Orchestrator = AccOrchestrator
    Linearizer = linearize

class DeviceCudaTarget(Target):
    Parizer = DeviceCudaizer
    DataManager = DeviceCudaDataManager
    Orchestrator = CudaOrchestrator
    Linearizer = cuda_linearize
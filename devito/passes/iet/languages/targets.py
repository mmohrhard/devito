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


class CTarget(Target):
    Parizer = SimdOmpizer
    DataManager = CDataManager


class OmpTarget(Target):
    Parizer = Ompizer
    DataManager = OmpDataManager


class DeviceOmpTarget(Target):
    Parizer = DeviceOmpizer
    DataManager = DeviceOmpDataManager
    Orchestrator = OmpOrchestrator


class DeviceAccTarget(Target):
    Parizer = DeviceAccizer
    DataManager = DeviceAccDataManager
    Orchestrator = AccOrchestrator

class DeviceCudaTarget(Target):
    Parizer = DeviceCudaizer
    DataManager = DeviceCudaDataManager
    Orchestrator = CudaOrchestrator
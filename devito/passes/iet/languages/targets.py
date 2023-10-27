from devito.passes.iet import linearize
from devito.passes.iet.languages.C import CDataManager
from devito.passes.iet.languages.openmp import (SimdOmpizer, Ompizer, DeviceOmpizer,
                                                OmpDataManager, DeviceOmpDataManager,
                                                OmpOrchestrator)
from devito.passes.iet.languages.openacc import (DeviceAccizer, DeviceAccDataManager,
                                                 AccOrchestrator)

__all__ = ['CTarget', 'OmpTarget', 'DeviceOmpTarget', 'DeviceAccTarget']


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

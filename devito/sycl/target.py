from devito.passes.iet.languages.targets import Target
from devito.sycl.passes.linearization import sycl_linearize


class DeviceSyclTarget(Target):
    from devito.sycl.passes import DeviceSyclizer, DeviceSyclDataManager, SyclOrchestrator

    Parizer = DeviceSyclizer
    DataManager = DeviceSyclDataManager
    Orchestrator = SyclOrchestrator
    Linearizer = sycl_linearize

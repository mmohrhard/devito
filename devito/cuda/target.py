from devito.passes.iet.languages.targets import Target


class DeviceCudaTarget(Target):
    from devito.cuda.passes import DeviceCudaizer, DeviceCudaDataManager, CudaOrchestrator

    Parizer = DeviceCudaizer
    DataManager = DeviceCudaDataManager
    Orchestrator = CudaOrchestrator

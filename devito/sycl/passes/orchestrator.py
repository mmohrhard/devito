from devito.passes.iet.orchestration import Orchestrator
from devito.sycl.lang import SyclBB


__all__ = ["SyclOrchestrator"]


class SyclOrchestrator(Orchestrator):
    lang = SyclBB

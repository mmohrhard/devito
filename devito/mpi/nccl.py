import ctypes
import sys

from devito.logger import error, info
from devito.types import Object

__all__ = ["NcclCommunicator", "NcclComm"]

from mpi4py import MPI

ncclUniqueId_t = ctypes.c_byte * 128
ncclCommunicator_t = ctypes.c_void_p


class NcclUniqueId(ctypes.Structure):
    def __init__(self, buf=None):
        super().__init__()
        if buf:
            self._opaque = buf

    _fields_ = [("_opaque", ncclUniqueId_t)]

    def __str__(self):
        return bytes(self._opaque).hex()


class NcclComm(Object):
    name = "nccl_comm"

    __rargs__ = ()

    def __init__(self, value=None):
        super().__init__(name=NcclComm.name, dtype=ctypes.c_void_p, value=value)

    def _arg_values(self, *args, **kwargs):
        grid = kwargs.get("grid", None)
        # Update `nccl_comm` based on object attached to `grid`
        if grid is not None:
            return grid.distributor._obj_nccl._arg_defaults()
        else:
            return self._arg_defaults()

    @property
    def _C_typename(self):
        return "ncclComm_t"


class NcclCommunicator:
    attempted_load = False
    nccl_lib = None

    @classmethod
    def get_nccl(cls):
        cls.attempted_load = True
        nccl_handle = "libnccl.so"
        try:
            from cuda.bindings.driver import CUresult, cuDeviceGetCount

            # We can't use NCCL if CUDA doesn't load properly
            ret, devcount = cuDeviceGetCount()
            if ret != CUresult.CUDA_SUCCESS:
                cls.nccl_lib = None
                return

            cls.nccl_lib = ctypes.CDLL(nccl_handle)
            cls.nccl_lib.ncclGetUniqueId.argtypes = [
                ctypes.POINTER(NcclUniqueId)
            ]
            cls.nccl_lib.ncclCommInitRank.argtypes = [
                ncclCommunicator_t,
                ctypes.c_int32,
                NcclUniqueId,
                ctypes.c_int32,
            ]

        except OSError:
            pass

        except ImportError:
            pass

    @classmethod
    def is_available(cls) -> bool:
        if not cls.attempted_load:
            cls.get_nccl()

        return cls.nccl_lib is not None

    @classmethod
    def should_use_for_mode(cls, mode: str | bool | None) -> bool:
        return isinstance(mode, str) and mode.startswith("nccl")

    def __init__(self, comm: MPI.Cartcomm):
        from cuda.bindings.driver import CUresult, cuDeviceGetCount
        from cuda.bindings.runtime import cudaSetDevice

        self._comm = comm
        info(
            f"attempting to set up a NCCL communicator over the top of MPI on rank {comm.rank}"
        )

        assert self.nccl_lib

        uid = NcclUniqueId()
        if comm.rank == 0:
            self.nccl_lib.ncclGetUniqueId(ctypes.byref(uid))

        # Broadcast the NCCL UID to all MPI ranks, then initialise the communicator
        comm.Bcast([uid, MPI.CHAR], 0)

        self._nccl_uid = uid

        global_rank = comm.rank

        split_comm = comm.Split_type(MPI.COMM_TYPE_SHARED, 0)
        local_rank = split_comm.rank

        ret, gpu_count = cuDeviceGetCount()

        if ret == CUresult.CUDA_SUCCESS and gpu_count > 0:
            info(
                f"binding MPI rank {global_rank} of {comm.size} to gpu {local_rank % gpu_count} of {gpu_count} on the local machine"
            )

            cudaSetDevice(local_rank % gpu_count)

            self._nccl_comm_t = ncclCommunicator_t(0)
            self.nccl_lib.ncclGroupStart()
            ret = self.nccl_lib.ncclCommInitRank(
                ctypes.byref(self._nccl_comm_t), comm.size, uid, comm.rank
            )
            self.nccl_lib.ncclGroupEnd()

            self._nccl_comm_obj = NcclComm(self._nccl_comm_t.value)

            if ret != 0:
                error(f"error setting up NCCL: {ret}")
                sys.exit(1)

            info("NCCL setup complete")

    def __del__(self):
        if self.nccl_lib:
            self.nccl_lib.ncclCommDestroy(self._nccl_comm_t)

    @property
    def comm_object(self) -> NcclComm:
        return self._nccl_comm_obj

from devito.logger import info, error
from devito.mpi.distributed import MPICommObject
from devito.types import Scalar, Object
import ctypes
import os

__all__ = ['NcclCommunicator', 'NcclComm']

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
    def __init__(self, name, value=None):
        super().__init__(name=name, dtype=ctypes.c_void_p, value=value)

    @property
    def _C_typename(self):
        return "ncclComm_t"
    

class NcclCommunicator:
    attempted_load = False

    @classmethod
    def get_nccl(cls):
        cls.attempted_load = True
        handle = 'libcudart.so'
        nccl_handle = '/d/sw/nvidia/hpc_sdk/22.9/Linux_x86_64/22.9/comm_libs/11.7/nccl/lib/libnccl.so.2.13.4'
        try:
            cls.cuda_lib = ctypes.CDLL(handle)

            c_devcount = ctypes.c_ulong(0)
            ret = cls.cuda_lib.cudaGetDeviceCount(ctypes.byref(c_devcount))
            if ret != 0 or c_devcount == 0:
                cls.cuda_lib = None

            cls.nccl_lib = ctypes.CDLL(nccl_handle)
            cls.nccl_lib.ncclGetUniqueId.argtypes = [ctypes.POINTER(NcclUniqueId)]
            cls.nccl_lib.ncclCommInitRank.argtypes = [ncclCommunicator_t, ctypes.c_int32, NcclUniqueId, ctypes.c_int32]

        except OSError:
            cls.cuda_lib = None
            cls.nccl_lib = None

    @classmethod
    def is_available(cls) -> bool:
        if cls.attempted_load == False:
            cls.get_nccl()

        return cls.cuda_lib is not None and cls.nccl_lib is not None

    def __init__(self, comm: MPI.Cartcomm):
        self._comm = comm
        info(f"attempting to set up a NCCL communicator over the top of MPI on rank {comm.rank}")

        uid = NcclUniqueId()
        if comm.rank == 0:
            self.nccl_lib.ncclGetUniqueId(ctypes.byref(uid))

        # Broadcast the NCCL UID to all MPI ranks, then initialise the communicator
        comm.Bcast([uid, MPI.CHAR], 0)
        
        self._nccl_uid = uid

        global_rank = comm.rank

        split_comm = comm.Split(MPI.COMM_TYPE_SHARED, 0)
        local_rank = split_comm.rank

        gpu_count = ctypes.c_long(0)
        ret = self.cuda_lib.cudaGetDeviceCount(ctypes.byref(gpu_count))
        if ret == 0 and gpu_count.value > 0:
            info(f"binding MPI rank {global_rank} of {comm.size} to gpu {local_rank % gpu_count.value} of {gpu_count.value} on the local machine")
            ret = self.cuda_lib.cudaSetDevice(local_rank % gpu_count.value)

            if ret != 0:
                error(f"error setting CUDA device: {ret}")
                os.exit(1)

            self._nccl_comm_t = ncclCommunicator_t(0)
            self.nccl_lib.ncclGroupStart()
            ret = self.nccl_lib.ncclCommInitRank(ctypes.byref(self._nccl_comm_t), comm.size, uid, comm.rank)
            self.nccl_lib.ncclGroupEnd()

            self._nccl_comm_obj = NcclComm("nccl_comm", self._nccl_comm_t.value)

            if ret != 0:
                error(f"error setting up NCCL: {ret}")
                os.exit(1)

            info("NCCL setup complete")        

    def __del__(self):
        ret = self.nccl_lib.ncclCommDestroy(self._nccl_comm_t)

    @property
    def comm_object(self) -> NcclComm:
        return self._nccl_comm_obj


NCCL_WORLD = None

from devito.arch.compiler import Compiler
import os

__all__ = ["SyclCompiler"]


class SyclCompiler(Compiler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.cflags.remove("-std=c99")
        self.cflags += ["-fsycl", "-std=c++17"]

        self.src_ext = "cpp"

    def __lookup_cmds__(self):
        # OneAPI Base Kit comes with dpcpp/icpx, both are clang++,
        # and icx, which is clang
        self.CC = "icx"
        self.CXX = "icpx"
        self.MPICC = "mpic++"
        self.MPICXX = "mpicxx"

    # only include dependencies that exist on disk - the SYCL compiler
    # generates and deletes some temporary files and they're gone by the
    # time codepy tries to look at them
    def get_dependencies(self, source_files):
        result = set(
            dep for dep in super().get_dependencies(source_files) if os.path.exists(dep)
        )
        return result

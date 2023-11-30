from spack.package import *


class Korc(CMakePackage):
    """KORC - Kinetic Orbit Runaway Code"""

    homepage = "https://ornl-fusion.github.io/KORC"
    url = "https://github.com/ORNL-Fusion/KORC/archive/refs/tags/v0.0.1.tar.gz"

    maintainers("mbeidler3", "cianciosa")

    license("UNKNOWN")

    version("0.0.1", sha256="b2eb25f42ed428250ad22e7462c628cb630e17ba4c08f38745161b0acffd9e2f")

    depends_on("gcc@13.1.0")
    depends_on("hdf5+fortran+mpi")

    def cmake_args(self):
        # FIXME: Add arguments other than
        # FIXME: CMAKE_INSTALL_PREFIX and CMAKE_BUILD_TYPE
        # FIXME: If not needed delete this function
        args = []
        return args


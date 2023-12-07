class Korc(CMakePackage):
    """KORC - Kinetic Orbit Runaway Code"""

    homepage = "https://ornl-fusion.github.io/KORC"
    url = "https://github.com/ORNL-Fusion/KORC/archive/refs/tags/v0.0.1.tar.gz"

    maintainers("mbeidler3", "cianciosa")

    license("UNKNOWN")

    version("0.0.1", sha256="b2eb25f42ed428250ad22e7462c628cb630e17ba4c08f38745161b0acffd9e2f")

    depends_on("gcc@13.1.0")
    depends_on("cmake %gcc@13.1.0")
    depends_on("hdf5+fortran+mpi %gcc@13.1.0")

    def cmake_args(self):
        args = [
            '-DUSE_PSPLINE=OFF'
            '-DUSE_FIO=OFF'
            '-DCMAKE_Fortran_FLAGS="-O3 -DHDF5_DOUBLE_PRESICION -fopenmp -malign-double -fconvert=\'big-endian\'"'
            '-DCMAKE_C_FLAGS="-O3 -fopenmp -malign-double"'
            '-DCMAKE_CXX_FLAGS="-O3 -fopenmp -malign-double"'
            '-DCMAKE_Fortran_FLAGS_DEBUG="-g -ffpe-trap=invalid,zero,overflow -fbacktrace -Werror"'
            '-DCMAKE_C_FLAGS_DEBUG="-g -g3"'
            '-DCMAKE_CXX_FLAGS_DEBUG="-g -g3"'
        ]
        return args

    def install(self, spec, prefix):
        mkdirp(prefix.bin)
        install('../spack-build/korc', prefix.bin)

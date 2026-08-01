#!/bin/bash

BUILD_TYPE=Debug

rm -f CMakeCache.txt
rm -rf CMakeFiles

rm -rf ./build_nvhpc_gpu && mkdir $_
cd build_nvhpc_gpu

  cmake \
    -DCMAKE_BUILD_TYPE:String=$BUILD_TYPE \
    -DUSE_OMP=OFF \
    -DUSE_FIO=OFF \
    -DUSE_ACC=ON \
    -DUSE_PSPLINE=ON \
    -DHDF5_PRECISION="DOUBLE" \
    -DCMAKE_Fortran_FLAGS="-acc=gpu -gpu=deepcopy,cc80 -Mfree -fPIC -c++libs -Mpreprocess -DDOUBLE_PRECISION" \
    -DCMAKE_Fortran_FLAGS_DEBUG="-O1" \
    -DCMAKE_C_FLAGS="-mp -DDOUBLE_PRECISION" \
    -DCMAKE_C_FLAGS_DEBUG="-g -G -traceback" \
    -DCMAKE_CXX_FLAGS="-std=c++11 -mp -DDOUBLE_PRECISION" \
    -DCMAKE_CXX_FLAGS_DEBUG="-g -traceback" \
    ..

make -j VERBOSE=1

#    -DCMAKE_Fortran_FLAGS_DEBUG="-O1 -g -Minstrument -traceback -lnvhpcwrapnvtx -gpu=debug,lineinfo" \
#using build_type=Regular or build_type=Debug with greater than -O1 leads to incorrect orbit calculations
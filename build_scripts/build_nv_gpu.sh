#!/bin/bash

BUILD_TYPE=Debug

rm -f CMakeCache.txt
rm -rf CMakeFiles

rm -rf ./build && mkdir $_
cd build

  cmake \
    -DCMAKE_BUILD_TYPE:String=$BUILD_TYPE \
    -DUSE_OMP=OFF \
    -DUSE_ACC=ON \
    -DUSE_PSPLINE=ON \
    -DCMAKE_Fortran_FLAGS="-acc=gpu -gpu=deepcopy,debug -Mfree -fPIC -c++libs -Mpreprocess -DHDF5_DOUBLE_PRESICION" \
    -DCMAKE_Fortran_FLAGS_DEBUG="-g -Minstrument -traceback -lnvhpcwrapnvtx" \
    -DCMAKE_C_FLAGS="-mp" \
    -DCMAKE_C_FLAGS_DEBUG="-g -G -traceback" \
    -DCMAKE_CXX_FLAGS="-std=c++11 -mp" \
    -DCMAKE_CXX_FLAGS_DEBUG="-g -traceback" \
    ..

make -j VERBOSE=1

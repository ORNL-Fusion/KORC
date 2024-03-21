#!/bin/bash

BUILD_TYPE=Debug

rm -f CMakeCache.txt
rm -rf CMakeFiles

rm -rf ./build && mkdir $_
cd build

cmake \
  -DCMAKE_BUILD_TYPE:String=$BUILD_TYPE \
  -DUSE_OMP=ON \
  -DUSE_PSPLINE=ON \
  -DCMAKE_Fortran_FLAGS="-Mfree -Mpreprocess -fPIC -O3 -mp -Mvect=simd:256 -c++libs -Mbyteswapio -DHDF5_DOUBLE_PRESICION" \
  -DCMAKE_Fortran_FLAGS_DEBUG="-g -Minstrument -traceback -lnvhpcwrapnvtx" \
  -DCMAKE_C_FLAGS="-O3 -mp -Mvect=simd:256" \
  -DCMAKE_C_FLAGS_DEBUG="-g -traceback" \
  -DCMAKE_CXX_FLAGS="-O3 -std=c++11 -mp -Mvect=simd:256" \
  -DCMAKE_CXX_FLAGS_DEBUG="-g -traceback" \
  ..

make -j VERBOSE=1

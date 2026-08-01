#!/bin/bash

#BUILD_TYPE=Release
BUILD_TYPE=Debug

rm -f CMakeCache.txt
rm -rf CMakeFiles

rm -rf ./build && mkdir $_
cd build

cmake \
  -DCMAKE_BUILD_TYPE:String=$BUILD_TYPE \
  -DUSE_FIO=OFF \
  -DUSE_PSPLINE=ON \
  -DUSE_ACC=ON \
  -DHDF5_PRECISION="DOUBLE" \
  -DCMAKE_Fortran_FLAGS="-acc=gpu -gpu=deepcopy,cc89 -Mfree -fPIC -c++libs -Mpreprocess -DDOUBLE_PRECISION" \
  -DCMAKE_C_FLAGS="-mp -DDOUBLE_PRECISION" \
  -DCMAKE_CXX_FLAGS="-std=c++11 -mp -DDOUBLE_PRECISION" \
  -DCMAKE_Fortran_FLAGS_DEBUG='-g -gpu=debug -O1 -Minfo=accel' \
  ..

make -j VERBOSE=1

#if [ $? -eq 0 ]
#  ctest --output-on-failure
#fi

#-DCMAKE_Fortran_FLAGS_DEBUG='-g -gpu=debug -O0 -Minfo=accel' \
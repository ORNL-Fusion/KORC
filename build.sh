#!/bin/bash

BUILD_TYPE=Debug

rm -f CMakeCache.txt
rm -rf CMakeFiles

rm -rf ./build && mkdir $_
cd build

cmake -DCMAKE_BUILD_TYPE:String=$BUILD_TYPE \
      -DUSE_PSPLINE=OFF \
      -DUSE_FIO=OFF \
      -DKORC_TEST=OFF \
      -DCMAKE_Fortran_FLAGS="-DHDF5_DOUBLE_PRESICION -malign-double -fconvert='big-endian'" \
      -DCMAKE_C_FLAGS="-malign-double"  \
      -DCMAKE_CXX_FLAGS="-malign-double" \
      -DCMAKE_Fortran_FLAGS_DEBUG="-ffpe-trap=zero,overflow -fbacktrace -Werror" \
      -DCMAKE_C_FLAGS_DEBUG="-g3" \
      -DCMAKE_CXX_FLAGS_DEBUG="-g3" ../

make -j VERBOSE=1

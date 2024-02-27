#!/bin/bash

BUILD_TYPE=Debug

rm -f CMakeCache.txt
rm -rf CMakeFiles

rm -rf ./build && mkdir $_


cmake -DCMAKE_BUILD_TYPE:String=$BUILD_TYPE \
      -DUSE_OMP=OFF \
      -DUSE_PSPLINE=OFF \
      -DUSE_FIO=OFF \
      -DCORI_DIR=OFF \
      -DKORC_TEST=OFF \
      -DCMAKE_Fortran_FLAGS="-O3 -DHDF5_DOUBLE_PRESICION -fopenmp -malign-double -fconvert='big-endian'" \
      -DCMAKE_C_FLAGS="-O3 -fopenmp -malign-double"  \
      -DCMAKE_CXX_FLAGS="-O3 -fopenmp -malign-double" \
      -DCMAKE_Fortran_FLAGS_DEBUG="-g -ffpe-trap=zero,overflow -fbacktrace -Werror" \
      -DCMAKE_C_FLAGS_DEBUG="-g -g3" \
      -DCMAKE_CXX_FLAGS_DEBUG="-g -g3"

make VERBOSE=1

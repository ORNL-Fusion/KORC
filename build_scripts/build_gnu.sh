#!/bin/bash

BUILD_TYPE=Debug

rm -f CMakeCache.txt
rm -rf CMakeFiles

rm -rf ./build && mkdir $_
cd build

export FC=/usr/bin/gfortran-11

cmake -DCMAKE_BUILD_TYPE:String=$BUILD_TYPE \
      -DUSE_PSPLINE=ON \
      -DUSE_FIO=OFF \
      -DKORC_TEST=OFF \
      -DCMAKE_Fortran_FLAGS="-malign-double -fconvert='big-endian'" \
      -DCMAKE_C_FLAGS="-malign-double"  \
      -DCMAKE_CXX_FLAGS="-malign-double" \
      -DCMAKE_Fortran_FLAGS_DEBUG="-g3 -ffpe-trap='zero,overflow' -ffpe-summary='all'  -fbacktrace" \
      -DCMAKE_C_FLAGS_DEBUG="-g3" \
      -DCMAKE_CXX_FLAGS_DEBUG="-g3" \
      -DFIO_LIBRARIES=/home/21b/KORC/FIO/install/lib/libfusionio.a  \
      -DM3DC1_LIBRARIES=/home/21b/KORC/FIO/install/lib/libm3dc1.a	\
      -DFIO_INCLUDE_PATH=/home/21b/KORC/FIO/install/include  \
      -DM3DC1_LIBRARIES="/home/21b/KORC/FIO/install/lib/libm3dc1.a;-L/home/21b/KORC/korc_gnu/.spack-env/view/lib;-lhdf5_fortran;-lhdf5;-lz;-ldl;-lm" \
      ../

make -j VERBOSE=1

#ctest --output-on-failure 

#testexit=$?

#if [ $testexit -eq 0 ]
#then
#  exit 0
#else
#  exit 1
#fi
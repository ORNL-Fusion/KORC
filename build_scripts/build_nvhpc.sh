#!/bin/bash

BUILD_TYPE=Debug

rm -f CMakeCache.txt
rm -rf CMakeFiles

rm -rf ./build && mkdir $_
cd build

cmake \
     -DUSE_PSPLINE=ON \
      -DUSE_ACC=OFF \
      -DUSE_OMP=OFF \
      -DUSE_FIO=OFF \
      -DKORC_TEST=OFF \
      -DCMAKE_Fortran_FLAGS="-DHDF5_DOUBLE_PRESICION -acc=multicore -c++libs" \
      -DCMAKE_C_FLAGS="-mp"   \
      -DCMAKE_CXX_FLAGS="-std=c++11 -mp" \
      -DCMAKE_CXX_FLAGS_DEBUG="-g -traceback -Bdynamic" \
      -DCMAKE_C_FLAGS_DEBUG="-g -traceback -Bdynamic" \
      -DCMAKE_Fortran_FLAGS_DEBUG="-g -Minfo=all -Minstrument -traceback -lnvhpcwrapnvtx" \
      -DFIO_LIBRARIES=/home/21b/KORC/FIO/install/lib/libfusionio.a  \
      -DM3DC1_LIBRARIES=/home/21b/KORC/FIO/install/lib/libm3dc1.a	\
      -DFIO_INCLUDE_PATH=/home/21b/KORC/FIO/install/include  \
      -DM3DC1_LIBRARIES="/home/21b/KORC/FIO/install/lib/libm3dc1.a;-L/home/21b/KORC/korc_gnu/.spack-env/view/lib;-lhdf5_fortran;-lhdf5;-lz;-ldl;-lm" \
      ../

make -j VERBOSE=1 

#ctest --output-on-failure
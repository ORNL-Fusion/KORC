#!/bin/bash

which h5diff
which nvfortran
which cmake
module list

cd ./KORC
rm -f CMakeCache.txt
rm -rf CMakeFiles
rm -rf ./build_gpu && mkdir $_
cd ./build_gpu
cmake \
    -DCMAKE_BUILD_TYPE:String=Debug \
    -DUSE_PSPLINE=ON \
    -DUSE_ACC=ON \
    -DUSE_FIO=OFF \
    -DCMAKE_Fortran_FLAGS="-acc=gpu -gpu=deepcopy -c++libs" \
    -DCMAKE_CXX_FLAGS="-std=c++11 -mp" \
    -DCMAKE_C_FLAGS="-mp" \
  ../
make -j VERBOSE=1

#
#    -gpu=cc80 indicates that target gpu is Nvidia A100
#

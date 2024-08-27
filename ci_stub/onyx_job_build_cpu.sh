#!/bin/bash

cd ./KORC
rm -f CMakeCache.txt
rm -rf CMakeFiles
rm -rf ./build_cpu && mkdir $_
cd ./build_cpu
cmake -DCMAKE_BUILD_TYPE:String=Debug \
    -DUSE_PSPLINE=ON \
    -DUSE_FIO=OFF \
  ../
make -j VERBOSE=1


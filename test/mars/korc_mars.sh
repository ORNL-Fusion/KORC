#!/bin/bash

set -ex

#define input file
INPUT_FILE=$1/"input_file_D3D_191366_1762ms_MARS.korc"
#define output directory
OUT_DIR="mars_test"

#check that output directory doesn't exist so bash doesn't complain
if [ ! -d $OUT_DIR ]; then
  mkdir -p $OUT_DIR
fi
if [ ! -f D3D_191366_1762ms_MARS.h5 ]; then
  ln -s $1/D3D_191366_1762ms_MARS.h5 D3D_191366_1762ms_MARS.h5
fi

./xkorc $INPUT_FILE $OUT_DIR/

h5diff -r -d 0.007 $OUT_DIR/file_0.h5 $1/file_0_new_random.h5


#!/bin/bash

#define input file
INPUT_FILE=$1/"input_file_egyro.korc"
#define output directory
OUT_DIR="egyro_test"

#check that output directory doesn't exist so bash doesn't complain
mkdir -p $OUT_DIR

#assumes binary directory ../KORC/build/bin was added to path
./xkorc $INPUT_FILE $OUT_DIR/

if [ $? == 0 ]; then
   echo "egyro test succeeded"
else
   echo "egyro test failed"
   exit $?
fi

h5diff -r -d 0.000001 $OUT_DIR/file_0.h5 $1/file_0.h5


#!/bin/bash

#define input file
INPUT_FILE="input_file_analytic_FO.korc"
#define output directory
OUT_DIR="OUT_TEST"

#check that output directory doesn't exist so bash doesn't complain
if [ ! -d "$OUT_DIR" ]; then
    mkdir $OUT_DIR
fi

#run KORC with input file outputting in desired directory
echo Running KORC with input file $INPUT_FILE and outputting to $OUT_DIR.

#assumes binary directory ../KORC/build/bin was added to path
xkorc $INPUT_FILE $OUT_DIR/ &

# pause before reading output file
sleep .25




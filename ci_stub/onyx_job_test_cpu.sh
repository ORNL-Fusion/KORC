#!/bin/bash

env | grep JOB_COUNT

JOB_COUNT=1

cd ./KORC/build_cpu && ctest --output-on-failure


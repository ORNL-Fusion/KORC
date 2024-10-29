#!/bin/bash

env | grep JOB_COUNT

JOB_COUNT=1

cd ./KORC/build_gpu && ctest -j ${JOB_COUNT:-1} --output-on-failure

for j in "mars_test egyro_test"; do
  for i in $(ls -1 ./${j}/ | sort -k1.6n); do

    echo "."
    echo "."
    echo "."
    echo "."
    echo "."
    echo "."
    echo "onyx_job_test_gpu >> [Info] << ${j}/${i}/output.korc"
    echo "."
    echo "."
    echo "."
    echo "."
    echo "."
    echo "."
    cat ./${j}/${i}/output.korc
  done
done

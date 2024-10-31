#!/bin/bash

env | grep JOB_COUNT

JOB_COUNT=1

pushd ./KORC/build_gpu >/dev/null

ctest -j --output-on-failure

pushd ./bin >/dev/null

for j in "mars_test" "egyro_test"; do
  for i in $(ls -1 ./${j}/ | sort -k1.6n); do
    for x in {1..20}; do echo "."; done
    echo "onyx_job_test_gpu >> [Info] << ${j}/${i}/output.korc"
    for x in {1..5}; do echo "."; done
    cat ./${j}/${i}/output.korc
  done
done


popd >/dev/null
popd >/dev/null

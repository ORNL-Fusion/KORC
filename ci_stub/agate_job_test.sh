#!/bin/bash

pushd ./KORC/build >/dev/null

ctest --output-on-failure
testexit=$?

pushd ./bin >/dev/null

for j in "mars_test" "egyro_test"; do
  for i in $(ls -1 ./${j}/ | sort -k1.6n); do
    for x in {1..20}; do echo "."; done
    echo "agate_job_test >> [Info] << ${j}/${i}/output.korc"
    for x in {1..5}; do echo "."; done
    cat ./${j}/${i}/output.korc
  done
done

popd >/dev/null
popd >/dev/null

if [ $testexit -eq 0 ]
then
  exit 0
else
  exit 1
fi
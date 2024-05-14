### Continuous Integration

The `main` branch of the Kinetic Orbit Runaway electrons Code (KORC) software repository should successfully build and
run on the computers the project is most interested in linux distributions
[supported at ORNL](https://ornl.servicenowservices.com/kb?id=kb_article_view&sysparm_article=KB0100342). At the time of
writing, they are RHEL 9 and Ubuntu 22.04 LTS. We are interested in the [Spack package manager](https://spack.io/), which
we use on [Perlmutter system](https://docs.nersc.gov/systems/perlmutter/) at NERSC. While we're happy to accept pull 
requests to address issues found on MacOS and Windows Subsystem for Linux (WSL), this is not a current priority of the 
project.

We are only testing on recent AMD and Intel systems using the x86_64 instruction set. Please contact us if your 
workstation/cluster uses a different architecture and would like it added to our test suite.

CPU-only builds are run on [GitHub-hosted runners](https://docs.github.com/en/actions/using-github-hosted-runners/about-github-hosted-runners/about-github-hosted-runners#supported-runners-and-hardware-resources). 

GPU builds are run on `milano0`, part of the
[ORNL ACSR Experimental Computing Laboratory](https://docs.excl.ornl.gov/system-overview). This machine has to NVIDIA A100
GPUS. Currents builds utilize the NVHPC compiler with OpenACC.

### Unit Testing

Unit testing is performed with the [FORTRAN Unit Test Framework (FRUIT)](https://sourceforge.net/projects/fortranxunit/). 
[xtest target](https://github.com/ORNL-Fusion/KORC/blob/18225f67b6c763fbfdf822f9a9a92f16cba39e00/test/CMakeLists.txt#L4)
These are being added retroactively primarily focused on new work.

### End-to-end Testing

Layer above the unit tests, are tests using the korc binary. All files go into a new directory inside [KORC/test](https://github.com/ORNL-Fusion/KORC/tree/18225f67b6c763fbfdf822f9a9a92f16cba39e00/test).  

ctest then runs a shell script, such as the one from egyo_test below.

```
#!/bin/bash

set -ex

#define input file
INPUT_FILE=$1/"input_file_egyro.korc"
#define output directory
OUT_DIR="egyro_test"

#check that output directory doesn't exist so bash doesn't complain
mkdir -p $OUT_DIR

#assumes binary directory ../KORC/build/bin was added to path
./xkorc $INPUT_FILE $OUT_DIR/

h5diff -r -d 0.000001 $OUT_DIR/file_0.h5 $1/file_0.h5
```



### Further Reading

https://www.oreilly.com/library/view/working-effectively-with/0131177052/
https://martinfowler.com/books/refactoring.html

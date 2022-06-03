
# KORC


 the Kinetic Orbit Runaway electrons Code (KORC) follows relativistic electrons in general electric and magnetic fields under the full Lorentz force, collisions, and radiation losses.


## Quick-Start

1. Find your platform under the `.cmfkit/platform` directory
```bash
[user@localhost korc]$ tree .cmfkit/platform/

.cmfkit/platform/
├── github-vm-macos-latest
│   ├── match
│   └── setup
├── github-vm-ubuntu-latest
│   ├── match
│   └── setup
├── github-vm-windows-latest
│   ├── match
│   └── setup
├── metal-macos-12
│   ├── match
│   └── setup
├── readybake
│   └── Dockerfile
├── readybake-alpine
│   └── Dockerfile
├── readybake-centos
│   └── Dockerfile
├── readybake-debian
│   └── Dockerfile
└── readybake-ubuntu
    └── Dockerfile

```

2. Create a `cmfkit_vars` file at the base directory of the project

```bash
#!/bin/bash

CMFKIT_PLATFORM='metal-macos-12'
```

1. configure cmfkit using environment vars
1. run the launcher script


## Required Dependencies

- GFortran compiler
- CMake >= 3.13
- MPI



differentiate by:
- compute_env:   [ nersc github docker ]
- os_type:       [ ubuntu macos ]
- mpi_wrapper:   [ openmpi mpich intel ]
-

Set some environment variables and run the launcher script like this:

```

CMFKIT_SCENARIO='egyro' \
CMFKIT_PLATFORM='readybake' \
CMFKIT_READYBAKE_IMAGE_SEED='ubuntu:22.04' \
CMFKIT_READYBAKE_INTERACTIVE='true' \
\
&& ./.cmfkit_base/launch

```

## More Info

http://ornl-fusion.github.io/KORC


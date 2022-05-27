
# KORC


 the Kinetic Orbit Runaway electrons Code (KORC) follows relativistic electrons in general electric and magnetic fields under the full Lorentz force, collisions, and radiation losses.


 ## Required Dependencies

- GFortran compiler
- CMake >= 3.13
- MPI


## Quick-Start Install

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

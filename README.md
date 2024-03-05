
# KORC
[![CI-CPU](https://github.com/ORNL-Fusion/KORC/actions/workflows/spack.yml/badge.svg)](https://github.com/ORNL-Fusion/KORC/actions/workflows/spack.yml)

 the Kinetic Orbit Runaway electrons Code (KORC) follows relativistic electrons in general electric and magnetic fields under the full Lorentz force, collisions, and radiation losses.


## Quick-Start

Spack provides an easy way to setup the gfortran compiler toolchain and other dependencies without the headache of managing them through the OS.  Here's how to build KORC on a Linux workstation using Spack:

1. Grab the spack software using git, and source the `setup-env.sh` script from within the spack installation.

```bash
[user@localhost KORC]$ git clone --depth=100 --branch=releases/v0.21 https://github.com/spack/spack.git ~/spack
[user@localhost KORC]$ . ~/spack/share/spack/setup-env.sh
```

2. Create and activate a new spack environment

```bash
[KORC] [user@localhost KORC]$ spack env create -d .
[KORC] [user@localhost KORC]$ spack env activate -p -d .
```


3. Configure spack and bootstrap the compiler

```bash
[KORC] [user@localhost KORC]$ spack config add "config:install_tree:padded_length:128"

[KORC] [user@localhost KORC]$ spack install --no-cache --add gcc@13.1.0
[KORC] [user@localhost KORC]$ spack load gcc@13.1.0
[KORC] [user@localhost KORC]$ spack compiler find

[KORC] [user@localhost KORC]$ spack install --no-cache --add cmake %gcc@13.1.0
[KORC] [user@localhost KORC]$ spack load cmake %gcc@13.1.0

[KORC] [user@localhost KORC]$ spack concretize -f
[KORC] [user@localhost KORC]$ spack install --no-cache --add hdf5+fortran+mpi %gcc@13.1.0
[KORC] [user@localhost KORC]$ spack load hdf5+fortran+mpi %gcc@13.1.0
```


4. Verify gfortan, cmake and Build KORC

```bash
[KORC] [user@localhost KORC]$ which cmake
[KORC] [user@localhost KORC]$ which gfortran
[KORC] [user@localhost KORC]$ gfortran --version

[KORC] [user@localhost KORC]$ ./build.sh
```

5. Exit the spack environment

```bash
[KORC] [user@localhost KORC]$ spack env deactivate
```


## More Info

http://ornl-fusion.github.io/KORC



FROM ubuntu:22.04

RUN set -x \
  && apt-get update \
  && DEBIAN_FRONTEND=noninteractive \
     apt-get install -y --no-install-recommends tzdata \
  && apt-get install -y \
                build-essential \
                gfortran \
                libnetcdf-dev \
                libnetcdff-dev \
                libopenmpi-dev \
                openmpi-common \
                openmpi-bin \
                cmake

#ENV FC=gfortran
#ENV CXX=g++
#ENV NETCDF_FLAGS=" -D_NETCDF"
#ENV L_NETCDF=" -L/usr/lib/x86_64-linux-gnu/lib -lnetcdff"
#ENV I_NETCDF=" -I/usr/include"


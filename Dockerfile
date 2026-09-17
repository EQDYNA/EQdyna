# EQdyna container image.
#
# Built and pushed by .github/workflows/publish.yml on every `vX.Y.Z` tag, as
# ghcr.io/eqdyna/eqdyna:<tag> and :latest. There is no other way this image is
# produced -- see Docker.guide.md for why that matters (it replaces a
# `docker commit` workflow that could not be rebuilt from source).
#
# Dependencies below are transcribed from install-eqdyna.sh's ubuntu branch
# (install-eqdyna.sh:68-70,91), not re-derived: this is MPICH, not OpenMPI,
# matching what install-eqdyna.sh actually invokes on ubuntu.
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# gfortran explicit, not implicit: mpich's package only depends on
# libgfortran5 (the runtime .so), not the gfortran compiler that `mpif90`
# shells out to for every build.
# NO --no-install-recommends: `mpich` carries the compiler wrappers but
# leaves mpif.h to a RECOMMENDED dev package, so the flag silently drops it.
# This file's job is to reproduce install-eqdyna.sh -m ubuntu, which uses a
# plain `apt-get install` -- do not use a stricter dependency resolution than
# the path being reproduced. Re-add the flag only alongside an explicitly
# enumerated dev-package list.
RUN apt-get update && \
    apt-get install -y \
        git vim make mpich gfortran \
        libnetcdf-dev libnetcdff-dev \
        python3 python3-pip \
        ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# python3-jax is deliberately NOT installed here: jax is only needed for the
# e2e sweep's python-jax backend column (a CI/dev-gate concern), not the
# Fortran production solver this image ships, and it would roughly double
# image size. CI (.github/workflows/test.yml) installs jax for that reason;
# this image does not.
#
# No --break-system-packages: that flag does not exist on ubuntu:22.04's
# stock pip3 (22.0.2 -- added in pip 23.0.1), and ubuntu:22.04 has no PEP 668
# externally-managed-environment marker on its system Python, so a plain
# pip3 install needs no flag at all.
RUN pip3 install numpy netCDF4 matplotlib xarray

COPY . /opt/eqdyna
WORKDIR /opt/eqdyna

RUN ./install-eqdyna.sh -m ubuntu

ENV EQDYNAROOT=/opt/eqdyna
ENV PATH="${EQDYNAROOT}/bin:${EQDYNAROOT}/scripts:${PATH}"

CMD ["bash"]

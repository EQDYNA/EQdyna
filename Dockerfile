# EQdyna container image.
#
# Built and pushed by .github/workflows/publish.yml on every `vX.Y.Z` tag, as
# ghcr.io/eqdyna/eqdyna:<tag> and :latest. See Docker.guide.md.
#
# Dependencies below are transcribed from install-eqdyna.sh's ubuntu branch
# (install-eqdyna.sh:68-70,91): MPICH, not OpenMPI, matching what
# install-eqdyna.sh actually invokes on ubuntu.
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# gfortran explicit: mpich's package depends only on libgfortran5 (the
# runtime .so), not the gfortran compiler that `mpif90` shells out to.
# NO --no-install-recommends: mpich leaves mpif.h to a RECOMMENDED dev
# package, so the flag would silently drop it.
RUN apt-get update && \
    apt-get install -y \
        git vim make mpich gfortran \
        libnetcdf-dev libnetcdff-dev \
        python3 python3-pip \
        ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# python3-jax is deliberately NOT installed here: it is only needed for the
# e2e sweep's python-jax backend column, not the Fortran production solver
# this image ships, and it would roughly double image size.
#
# No --break-system-packages: ubuntu:22.04's stock pip3 (22.0.2) predates
# that flag (pip 23.0.1+) and has no PEP 668 marker on its system Python.
RUN pip3 install numpy netCDF4 matplotlib xarray

COPY . /opt/eqdyna
WORKDIR /opt/eqdyna

RUN ./install-eqdyna.sh -m ubuntu

ENV EQDYNAROOT=/opt/eqdyna
ENV PATH="${EQDYNAROOT}/bin:${EQDYNAROOT}/scripts:${PATH}"

CMD ["bash"]

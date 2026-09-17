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
# libgfortran5 (the runtime .so), not the gfortran compiler `mpif90` shells
# out to for every build -- `make` failed with "gfortran: command not
# found" INSIDE this exact RUN step (found 2026-09-16, this image's first
# real build after the CI-gate fixes landed). install-eqdyna.sh's ubuntu
# branch had the identical gap and is fixed the same way, same commit --
# every prior environment that ran it (this dev box, GitHub Actions'
# ubuntu-22.04 runner) already had gfortran from something else, so a bare
# ubuntu:22.04 base image is what finally exposed it.
# NO --no-install-recommends. It was here, and it is what made this image
# diverge from the very script this file exists to reproduce. `mpich` carries
# the compiler wrappers but leaves mpif.h to a RECOMMENDED dev package, so the
# flag silently removed it and the build died at
#   errorCodes.f90:122: Error: Can't open included file 'mpif.h'
# one layer after the gfortran fix -- the second missing-package failure in a
# row from the same cause. install-eqdyna.sh:77-79 and test.yml both run a
# plain `apt-get install` and both build fine; this file's whole job is to
# reproduce install-eqdyna.sh -m ubuntu, so it must not use a STRICTER
# dependency resolution than the path it is reproducing. Re-add the flag only
# alongside an explicitly enumerated dev-package list, verified by a
# `gh workflow run publish.yml` dispatch, as its own change.
RUN apt-get update && \
    apt-get install -y \
        git vim make mpich gfortran \
        libnetcdf-dev libnetcdff-dev \
        python3 python3-pip \
        ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# python3-jax is deliberately NOT installed here. jax is only needed for the
# e2e sweep's python-jax backend column, which is a CI/dev-gate concern, not
# something the Fortran production solver (what this image ships) needs to
# run. jax is also a large dependency (pulls in its own numpy/scipy pins) that
# would roughly double image size for a backend most container users will
# never touch. CI (.github/workflows/test.yml) installs jax for that reason;
# this image does not. Revisit if a jax-GPU container variant is added later
# (a CUDA variant would need its own separate Dockerfile, not this one).
#
# No --break-system-packages: that flag does not exist on ubuntu:22.04's
# stock pip3 (22.0.2 -- the flag was added in pip 23.0.1) and its absence
# here broke the FIRST real build this Dockerfile ever went through
# (v5.8.3's tag-triggered publish.yml run, exit code 2, "no such option");
# nothing had ever actually built this image before that run, since the
# original publish workflow shipped with zero build verification. ubuntu:22.04
# has no PEP 668 externally-managed-environment marker on its system Python,
# so a plain pip3 install needs no flag at all.
RUN pip3 install numpy netCDF4 matplotlib xarray

COPY . /opt/eqdyna
WORKDIR /opt/eqdyna

RUN ./install-eqdyna.sh -m ubuntu

ENV EQDYNAROOT=/opt/eqdyna
ENV PATH="${EQDYNAROOT}/bin:${EQDYNAROOT}/scripts:${PATH}"

CMD ["bash"]

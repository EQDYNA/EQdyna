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

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git vim make mpich \
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
# (see the note in QUEUE_dockerfile.md about a CUDA variant / separate
# Dockerfile, not this one).
RUN pip3 install --break-system-packages numpy netCDF4 matplotlib xarray

COPY . /opt/eqdyna
WORKDIR /opt/eqdyna

RUN ./install-eqdyna.sh -m ubuntu

ENV EQDYNAROOT=/opt/eqdyna
ENV PATH="${EQDYNAROOT}/bin:${EQDYNAROOT}/scripts:${PATH}"

CMD ["bash"]

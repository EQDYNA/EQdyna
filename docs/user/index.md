# EQdyna User Guide

EQdyna is a parallel finite element code for simulating earthquake
spontaneous dynamic rupture, seismic wave propagation, and high-frequency
deterministic ground motions on geometrically complex fault systems. The
solver core is written in Fortran, with a from-scratch Python/JAX port of
the same physics for verification and GPU use, and Python utilities for
setting up and running cases.

This site covers installing EQdyna, setting up and running a case, the
parameter reference, output file formats, the SCEC benchmark suite EQdyna
is verified against, and measured performance.

## Where to go

* Getting started -- installing EQdyna on Ubuntu, an HPC cluster, or macOS, and running a first case.
* Running a case -- creating a case, configuring it, choosing an MPI rank count, and the Python/JAX backend.
* Parameters -- every configurable default, generated from the code itself.
* Output files -- what each output file contains and how its columns are defined.
* Benchmarks -- the SCEC verification suite and how each case is checked.
* Performance -- measured wall-clock and memory figures, with the hardware they were measured on.
* Troubleshooting -- exit codes and what a failed run is telling you.
* Citing -- how to cite EQdyna.

## Getting the code

EQdyna is released under the MIT License and hosted on GitHub. A Docker
image is also available; see the repository's Docker guide for the
container-based route if you would rather not build from source.

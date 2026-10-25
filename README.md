# EQdyna

EQdyna is a parallel finite element code for simulating earthquake
spontaneous dynamic rupture, seismic wave propagation, and high-frequency
deterministic ground motions on geometrically complex fault systems. The
solver core is written in Fortran90, with Python utilities for setting up
cases and a from-scratch Python/JAX port of the same physics for
verification and GPU use. EQdyna is also the dynamic-rupture engine inside
the multicycle earthquake simulator EQsimu
([Liu et al., 2020, GJI](https://www.researchgate.net/publication/346814142_EQsimu_a_3-D_finite_element_dynamic_earthquake_simulator_for_multicycle_dynamics_of_geometrically_complex_faults_governed_by_rate-_and_state-dependent_friction)).

Released under the MIT License (`LICENSE`).

## Requirements

* A Fortran compiler (gfortran or Intel Fortran), MPI (mpich or Intel MPI) and netCDF (libnetcdf, libnetcdff).
* Python 3 with numpy>=1.20, matplotlib, xarray and netCDF4.

On Ubuntu 22 this is one root step, then a Python environment of your own. The
environment keeps EQdyna's Python packages apart from the system's, which may be
built for an older numpy:

```
sudo apt-get install git make gfortran mpich libnetcdf-dev libnetcdff-dev python3 python3-pip   # needs root
python3 -m pip install --user virtualenv
python3 -m virtualenv ~/eqdyna-env
. ~/eqdyna-env/bin/activate
pip install numpy netCDF4 matplotlib xarray jax
```

For other platforms (macOS, TACC Lonestar6, Texas A&M Grace) and GPU
support, see [docs/user/getting-started.md](docs/user/getting-started.md).

## Install

```
git clone https://github.com/EQDYNA/EQdyna.git
cd EQdyna
./install-eqdyna.sh -m ubuntu      # -e macos on a Mac; -m ls6 on TACC Lonestar6
export EQDYNAROOT=$(pwd)
export PATH=$EQDYNAROOT/bin:$EQDYNAROOT/scripts:$PATH
export PYTHONPATH=$EQDYNAROOT/src/python
```

Add the three `export` lines and `. ~/eqdyna-env/bin/activate` to your `~/.bashrc`
(with the full path in place of `$(pwd)`) so every new shell finds EQdyna. If you
later edit `src/fortran/`, re-run `./install-eqdyna.sh`; the tools refuse to use a
binary built from older source. The quick start below is the install check: if it
produces the files it lists, the install works.

## Quick start

Run the TPV8 benchmark in a new directory:

```
create.newcase ~/runs/tpv8 test.tpv8
cd ~/runs/tpv8
./case.setup
bash run.sh
```

`bash run.sh` runs 5 s of simulated time on 4 MPI ranks, in about 15 seconds on a
workstation. It leaves these results in the case directory:

* `cRuptureDynamics.png`, a plot of rupture time, slip and peak slip rate on the fault.
* `frt.txt*`, the on-fault rupture time, slip and stress at every fault node.
* `faultst*.txt` and `body*.txt`, time series at on-fault and off-fault stations.
* `fault.dyna.r.nc`, the final fault state in netCDF.

To run another case, replace `test.tpv8` with one of the SCEC/USGS benchmark
cases listed in [docs/user/benchmarks.md](docs/user/benchmarks.md). To build
your own, copy the closest case and edit its `user_defined_params.py`. On an
HPC cluster, `./case.submit` submits the run as a batch job instead.

### Python solver

The same physics is also implemented in Python/JAX. It runs a case on one
process, so first set the decomposition to one rank:

```
create.newcase ~/runs/tpv8-jax test.tpv8
cd ~/runs/tpv8-jax
printf '\npar.nx = 1\npar.ny = 1\npar.nz = 1\n' >> user_defined_params.py
./case.setup
python3 -m eqdyna . --backend jax
```

It writes `frt.txt0` and the station files in about 25 seconds. The default device
is the CPU; add `--device cuda` to run on a GPU.

A slower, spec-resolution reproduction of the full SCEC benchmark suite is
opt-in and can take hours: `EQDYNA_FULL_LAUNCH=yes-hours python3
testsys/run.py e2e-full`. See
[docs/user/benchmarks.md](docs/user/benchmarks.md) for what it checks.

## Where to go next

* [docs/user/index.md](docs/user/index.md) -- the full user guide: running a case, parameters, output files, benchmarks and performance.
* [docs/user/troubleshooting.md](docs/user/troubleshooting.md) -- exit codes and what they mean.
* [Docker.guide.md](Docker.guide.md) -- running EQdyna from a container instead of building it.
* [docs/user/citing.md](docs/user/citing.md) -- how to cite EQdyna, and how to reach the maintainers.

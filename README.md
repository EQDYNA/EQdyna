# EQdyna

EQdyna is a parallel finite element code for simulating earthquake
spontaneous dynamic rupture, seismic wave propagation, and high-frequency
deterministic ground motions on geometrically complex fault systems. The
solver core is written in Fortran90, with Python utilities for setting up
cases and a from-scratch Python/JAX port of the same physics for
verification and GPU use. EQdyna is also the dynamic-rupture engine inside
the multicycle earthquake simulator EQsimu
([Liu et al., 2020, GJI](https://www.researchgate.net/publication/346814142_EQsimu_a_3-D_finite_element_dynamic_earthquake_simulator_for_multicycle_dynamics_of_geometrically_complex_faults_governed_by_rate-_and_state-dependent_friction)).

## Latest release

* 20260924 v5.17.0 release notes. See the [GitHub Release](https://github.com/EQDYNA/EQdyna/releases/tag/v5.17.0) for the full notes.
Past releases are archived in `pastReleaseNotes.md`.

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

`ubuntu.env.sh` installs the same system packages and must be run as root. For
NVIDIA GPUs, install `"jax[cuda12]"` instead of `jax`. On macOS, install
[Homebrew](https://brew.sh) plus `brew install gcc netcdf netcdf-fortran`; the
`./install-eqdyna.sh -e macos` below then adds mpich and the Python packages.

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

To run another case, replace `test.tpv8` with one of the case names below. To
build your own, copy the closest case and edit its `user_defined_params.py`. On
an HPC cluster, `./case.submit` submits the run as a batch job instead.

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

## Benchmarks

EQdyna is verified against SCEC/USGS benchmarks. A fast local check runs
every case at coarse resolution (`python3 testsys/run.py e2e`); a slower, spec-resolution
reproduction is opt-in and can take hours
(`EQDYNA_FULL_LAUNCH=yes-hours python3 testsys/run.py e2e-full`).

| case | physics | SCEC benchmark | details |
|---|---|---|---|
| test.tpv8 | strike-slip, slip-weakening | [TPV8](https://strike.scec.org/cvws/tpv89docs.html) | [README](case_input/test.tpv8/README.md) |
| test.tpv10 | dipping normal fault | [TPV10](https://strike.scec.org/cvws/tpv10_11docs.html) | [README](case_input/test.tpv10/README.md) |
| test.tpv104 | strike-slip, rate-and-state friction | [TPV104](https://strike.scec.org/cvws/tpv103_104docs.html) | [README](case_input/test.tpv104/README.md) |
| test.tpv1053d | rate-and-state + thermal pressurization | [TPV105-3D](https://strike.scec.org/cvws/tpv105_3D_docs.html) | [README](case_input/test.tpv1053d/README.md) |
| test.tpv29 | strike-slip, fractal rough fault | [TPV29](https://strike.scec.org/cvws/tpv29_30docs.html) | [README](case_input/test.tpv29/README.md) |
| test.tpv30 | rough fault + off-fault plasticity | [TPV30](https://strike.scec.org/cvws/tpv29_30docs.html) | [README](case_input/test.tpv30/README.md) |
| test.drv.a6 | fractal fault + off-fault plasticity | internal, no published spec | [README](case_input/test.drv.a6/README.md) |

See `docs/user/benchmarks.md` for the full case list and how verification
works, and `docs/user/performance.md` for wall-clock and memory figures.

## Where to go next

* `docs/user/benchmarks.md` -- the full benchmark list and verification approach
* `docs/user/performance.md` -- measured speed and memory figures
* `docs/user/troubleshooting.md` -- exit codes and what they mean
* `Docker.guide.md` -- running EQdyna from a container instead of building it
* Each `case_input/<case>/README.md` -- details specific to that benchmark

## Collaboration

We try our best to make EQdyna easy to use, but doing great science is our
priority. We welcome collaborations, comments, and suggestions. Please
reach out to Drs Benchun Duan (bduan@tamu.edu) and Dunyu Liu
(dliu@ig.utexas.edu).

## Citation and references

EQdyna is released under the MIT License (`LICENSE`). If you use it, please
cite the relevant paper(s) below:

1. [Duan and Oglesby (2006)](https://doi.org/10.1029/2005JB004138).
   Heterogeneous fault stresses from previous earthquakes and the effect on
   dynamics of parallel strike-slip faults, J. Geophys. Res., 111.
2. [Duan (2012)](https://doi.org/10.1029/2011JB009124). Dynamic rupture of
   the 2011 Mw 9.0 Tohoku-Oki earthquake: Roles of a possible subducting
   seamount, J. Geophys. Res., 117.
3. [Luo and Duan (2018)](https://doi.org/10.1029/2017JB015320). Dynamics of
   non-planar thrust faults governed by various friction laws, J. Geophys.
   Res. Solid Earth, 123.
4. [Liu and Duan (2018)](https://doi.org/10.1785/0120170374). Scenario
   Earthquake and Ground-Motion Simulations in North China: Effects of
   Heterogeneous Fault Stress and 3D Basin Structure, Bull. Seismol. Soc.
   Am., 108(4).

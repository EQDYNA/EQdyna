# Getting Started

## Requirements

* A Fortran compiler (gfortran or Intel Fortran)
* MPI (mpich or Intel MPI)
* netCDF (libnetcdf, libnetcdff)
* Python 3 with numpy, matplotlib, xarray, netCDF4

## Install

Clone the repository first:

```
git clone https://github.com/EQDYNA/EQdyna.git
cd EQdyna
```

Then install for your machine. `install-eqdyna.sh` builds the Fortran
solver and moves the executable to `bin/`; on Ubuntu and macOS it can also
install the system and Python dependencies listed above.

* **Ubuntu 22.04**: `bash ubuntu.env.sh` installs the dependencies through
  `apt-get` and `pip`, then `./install-eqdyna.sh -m ubuntu` builds.
* **TACC Lonestar6**: `./install-eqdyna.sh -m ls6` loads the cluster's
  `netcdf` module and builds; no separate dependency install is needed.
* **Texas A&M Grace**: `./install-eqdyna.sh -m grace` loads the cluster's
  netCDF module and builds.
* **macOS**: `./install-eqdyna.sh -e macos` installs `mpich` and `python`
  through Homebrew plus the required Python packages, then builds.

After building, set the environment variables every session needs:

```
export EQDYNAROOT=$(pwd)
export PATH=$EQDYNAROOT/bin:$EQDYNAROOT/scripts:$PATH
export PYTHONPATH=$EQDYNAROOT/src/python
```

Add those lines to your shell's startup file so a new shell always has EQdyna
on its path. `PYTHONPATH` is what lets `python3 -m eqdyna` find the Python
solver. The build stamps `bin/eqdyna` with a hash of the Fortran source; after
editing `src/fortran/`, re-run `./install-eqdyna.sh`, because the test tools
refuse a binary built from older source. Optional, for running the Python solver on a GPU:
`pip install "jax[cuda12]"`.

Check the install with the fast test tiers:

```
python3 testsys/run.py unit regression
```

This needs `pip install pytest` in the same environment. It takes about 4
minutes (246 s measured on 2026-09-25 on a shared 64-core workstation) and must
end with `SUCCESS unit` and `SUCCESS regression`.

## Running a first case

The quickest way to see EQdyna run end to end is the SCEC TPV8 benchmark, a
vertical strike-slip fault in a homogeneous elastic halfspace:

```
create.newcase mytpv8 test.tpv8
cd mytpv8
./case.setup
bash run.sh
```

`run.sh` runs the solver under MPI and then produces a rupture-dynamics
plot, `cRuptureDynamics.png`, in the same directory. See Running a Case for
what each step does and how to configure your own case.

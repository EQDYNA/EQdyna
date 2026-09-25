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

* A Fortran compiler (gfortran or Intel Fortran)
* MPI (mpich or Intel MPI)
* netCDF (libnetcdf, libnetcdff)
* Python 3 with numpy>=1.20, matplotlib, xarray, netCDF4

On Ubuntu 22, `bash ubuntu.env.sh` installs all of the above through
apt-get and pip.

Optional, for the Python solver on a GPU: `pip install "jax[cuda12]"`, then
run with `python3 -m eqdyna <case_dir> --backend jax`; verify with
`python3 testsys/run.py gpu`.

## Install

```
git clone https://github.com/EQDYNA/EQdyna.git
cd EQdyna
chmod 755 install-eqdyna.sh
./install-eqdyna.sh -m ubuntu   # ubuntu / ls6 / macos
export EQDYNAROOT=$(pwd)
PATH=$EQDYNAROOT/bin:$EQDYNAROOT/scripts:$PATH
python3 testsys/run.py unit regression   # a few seconds
```

For bash, add these two lines to `.bashrc` so every new shell has EQdyna on
its path:

```
export EQDYNAROOT=/path/to/EQdynaRootDirectory
PATH=$EQDYNAROOT/bin:$EQDYNAROOT/scripts:$PATH
```

The build stamps `bin/eqdyna` with a hash of the Fortran source it was built
from. After editing anything under `src/fortran/`, rebuild with
`./install-eqdyna.sh`; the test tools otherwise refuse with
`bin/eqdyna was built from different source ... rebuild with ./install-eqdyna.sh`.

## Quick start

Three steps run a new case:

```
create.newcase $caseDirectoryName $predefinedCompset
cd $caseDirectoryName
./case.setup
bash run.sh   # or ./case.submit to submit a batch job on an HPC cluster
```

`$predefinedCompset` is one of the pre-defined cases below, or your own copy
of one with `user_defined_params.py` edited to match your fault and loading:

* test.drv.a6, for deterministic ground motion with a fractal fault and plasticity
* test.tpv8, test.tpv10, test.tpv36, test.tpv37, test.tpv104, test.tpv1053d, test.tpv29, test.tpv30

`TPV<number>` follows the naming convention of the [SCEC/USGS Spontaneous Rupture Code Verification Project](https://strike.scec.org/cvws/).

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

# News in 2026
* 20260914 v5.4.0 release notes
  * New - python/eqdyna/standalone: a fully standalone Python EQdyna (zero Fortran involved), now covering all five gated cases (tpv8, tpv104, tpv1053d, tpv10, drv.a6) with a committed acceptance tier (`testsys/run.py accept`) — four at roundoff-level agreement, drv.a6 under a documented chaos-aware criterion for its rupture-arrest bistability.
  * New - Drucker-Prager viscoplasticity ported to the standalone solver (test.drv.a6).
  * New - optional GPU execution (JAX/CUDA), verified via `testsys/run.py gpu`.
  * Add - testsys/ parity, accept, gpu, perf, and scaling tiers.
  * Fix - PML boundary tests standardized to inclusive bounds; mesh-time checkPMLAlignment guard added.
  * Fix - retired dead fric slot 6 (unused surface-pore-pressure reads).
  * Refactor - root cleanup: netCDF4 replaces xarray in check.test.py; orphaned testNameListWhole.py removed; shared loading/kernels modules factored out of the Python port.
  * Fix - scripts/lib.py: lazy imageio import, guarded optional dependencies.
  * Full details: the v5.4.0 GitHub Release, git tag message, and pathway_forward.md.
  * For past release notes, please refer to pastReleaseNotes.md.

# Introduction to *```EQdyna```*

*```EQdyna```* is a parallel finite element software to simulate earthquake spontaneous dynamic rupture, seismic wave propagation and high frequency deterministic ground motions. It has a focus to simulate earthquakes on geometrically complex fault systems with other heterogeneities in mind. The core is written in FORTRAN90 with a set of python utilities to setup cases and allocate HPC resources.  <br/>

*```EQdyna```* is highly efficient and scalable due to its adoptions of built-in mesh generation, explicit time integration, under-integrated hexahedral elements, degenerated wedge elements for complex fault geometries, hourglass controls, perfectly matched layer and 3D parallelization with MPI. <br/> 

Other features include 
* 3D velocity structure to simulate basin effects.
* Frequency independent Q by coarsed-grained memory scheme for seismic attenuation. 
* Various frictional constitutions such as slip-weakening and various forms of rate- and state-friction.
* Normal stress evolutions by an additional state variable.
* Drucker-prager off-fault viscoplasticity.
* Dynamic relaxation for earthquake cycle applications, etc.

*```EQdyna```* has been extensively verified against benchmark problems from [SCEC/USGS Spontaneous Rupture Code Verification Project](https://strike.scec.org/cvws/).

*```EQdyna```* is also part of the fully dynamic earthquake cycle simulator *```EQsimu```* [(*Liu et al., 2020, GJI*)](https://www.researchgate.net/publication/346814142_EQsimu_a_3-D_finite_element_dynamic_earthquake_simulator_for_multicycle_dynamics_of_geometrically_complex_faults_governed_by_rate-_and_state-dependent_friction).

# Verification & Benchmarks

EQdyna is verified against SCEC/USGS Spontaneous Rupture Code Verification
benchmarks (https://strike.scec.org/cvws/). The fast tier gates every
commit at coarse resolution against frozen references
(`python3 testsys/run.py e2e`); the full tier reproduces each gated
benchmark at its official spec resolution and duration, report-only, on
16 ranks (`python3 testsys/run.py e2e-full`).

| case | physics | SCEC benchmark | fast tier (dx/cores/~time) | full-tier spec (dx/term) |
|---|---|---|---|---|
| [test.tpv8](case_input/test.tpv8/README.md) | strike-slip, slip-weakening | [TPV8](https://strike.scec.org/cvws/tpv89docs.html) | 500m/4/~48s | 100m/15s |
| [test.tpv10](case_input/test.tpv10/README.md) | dipping normal fault | [TPV10](https://strike.scec.org/cvws/tpv10_11docs.html) | 500m/4/~48s | 100m/15s |
| [test.tpv104](case_input/test.tpv104/README.md) | strike-slip, rate-and-state | [TPV104](https://strike.scec.org/cvws/tpv103_104docs.html) | 500m/4/~48s | 50m/12s |
| [test.tpv1053d](case_input/test.tpv1053d/README.md) | RSF + thermal pressurization | [TPV105-3D](https://strike.scec.org/cvws/tpv105_3D_docs.html) | 500m/4/~48s | excluded (no spec dx) |
| [test.drv.a6](case_input/test.drv.a6/README.md) | fractal fault + plasticity | internal | 500m/4/~48s | excluded (no published spec) |
| [test.meng2023a](case_input/test.meng2023a/README.md) | layered velocity structure | internal | 400m/4/~48s | excluded (no published spec) |
| [test.meng2023cb](case_input/test.meng2023cb/README.md) | layered velocity, multi-patch | internal | 400m/4/~48s | excluded (no published spec) |

# Environment
*Optional (Python solver on GPU)*: `pip install "jax[cuda12]"` — the standalone Python solver (`python -m eqdyna.standalone`) then runs on NVIDIA GPUs; verify with `python3 testsys/run.py gpu`.

*```EQdyna```* requires <br/>
  - FORTRAN compiler (gfortran/intel FORTRAN)
  - MPI (mpich/intel MPI)
  - netCDF (libnetcdf libnetcdff)

Pre-staging and post-processing require Python packages <br/>
  - Python3
  - numpy>=1.20 (as of 20240201, matplotlib requires numpy>=1.20)
  - matplotlib
  - xarray
  - netCDF4
```
bash ubuntu.env.sh
```
will install the required packages through apt-get and pip on Ubuntu 22. <br/>

# Installation
```
git clone https://github.com/EQDYNA/EQdyna.git
cd EQdyna
chmod 755 install-eqdyna.sh
./install-eqdyna.sh -m ubuntu # ubuntu/ls6/macos
export EQDYNAROOT=$(pwd)
PATH=$EQDYNAROOT/bin:$EQDYNAROOT/scripts:$PATH
python3 testAll.py # quick testing of multiple examples with 4 cores; take ~180 s.
```
For bash, please insert the following lines in .bashrc
```
export EQDYNAROOT=/path/to/EQdynaRootDirectory
PATH=$EQDYNAROOT/bin:$EQDYNAROOT/scripts:$PATH
```

# Quick Start Guide
Three steps are needed to run a new case. <br/>
```
create.newcase $caseDirectoryName $predefinedCompset
cd $caseDirectoryName
./case.setup
bash run.sh # or ./case.submit to submit batch job on LS6.
```
Replace $caseDirectoryName with the directory name you want to create. <br/>
Replace $predefinedCompset with one of the following supported compsets <br/>
* test.drv.a6, for determinisitc ground motion with fractal fault and plasticity
* test.tpv8
* test.tpv10
* test.tpv36
* test.tpv37
* test.tpv104
* test.tpv1053d <br/>
[TPV+number is the naming convention of [SCEC/USGS Spontaneous Rupture Code Verification Excercise](https://strike.scec.org/cvws/).] <br/>

For a customized case, please choose the most relevant predefined compset and modify ```user_defined_param.py``` accordingly. <br/>

# Benchmark computational performance and resource
* TPV36 & 37: 4.7 hours for 50 m resolution using 512 CPUs on Lonestar6 at TACC. <br/>
* TPV104: 0.4 hours for 15-sec simulation (1875 time steps) using 40 CPUs on Lonestar6.  <br/>

# Collaboration
We try our best to make *```EQdyna```* easy to use but doing great science is our priority. We welcome collaborations, comments and suggestions. Please reach out to Drs Benchun Duan (bduan@tamu.edu) and Dunyu Liu (dliu@ig.utexas.edu).

# Reference
1. [Duan and Oglesby (2006)]( https://doi.org/10.1029/2005JB004138). Heterogeneous fault stresses from previous earthquakes and the effect on dynamics of parallel strike-slip faults, J. Geophys. Res., 111.
2. [Duan (2012)](https://doi.org/10.1029/2011JB009124). Dynamic rupture of the 2011 Mw 9.0 Tohoku-Oki earthquake: Roles of a possible subducting seamount, J. Geophys. Res., 117.
3. [Luo and Duan (2018)](https://doi.org/10.1029/2017JB015320). Dynamics of non-planar thrust faults governed by various friction laws, J. Geophys. Res. Solid Earth, 123.
4. [Liu and Duan (2018)](https://doi.org/10.1785/0120170374). Scenario Earthquake and Ground‐Motion Simulations in North China: Effects of Heterogeneous Fault Stress and 3D Basin Structure, Bull. Seismol. Soc. Am., 108(4).

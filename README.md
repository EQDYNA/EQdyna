# News in 2026
* 20260909 v5.3.5 release notes
  * Add - `testsys/`, a tiered test system (unit/regression/e2e) with a single entry point `testsys/run.py [unit|regression|e2e|all]`; `testAll.py` is now a thin wrapper delegating to `testsys/e2e/run_e2e.py` (PROJECT_RULES.md rule 3). Gate: `python3 testsys/run.py unit` (33 tests, pure Python, no MPI/Fortran, <1s) then `regression` (2 guards, one per past incident, rule 10) then `e2e` (full create.newcase -> case.setup -> mpirun -> plotRuptureDynamics -> check.test.py pipeline against test.reference.results/, rule 7).
  * Fix - `scripts/lib.py`: `B2`/`B3`'s negative-`y` guard called `sys.exit()`, but the file did `from sys import *`, which does not bind the name `sys` -- the guard raised `NameError` instead of exiting. Changed to `import sys`. Covered by new unit tests `test_B2_negative_y_exits_cleanly_not_with_nameerror` / `test_B3_...` in `testsys/unit/test_lib.py`.
  * Fix - `check.test.py`'s `compare_nc_files` set `metadata_equal = f1.identical(f2)`, which requires bit-exact data in addition to matching attrs -- conflating the function's one calibrated 1e-3 threshold (rule 5) with an uncalibrated bit-exact check. This turned ordinary parallel-MPI floating-point non-determinism into a false `FAIL metadata` (observed: `test.tpv10/fault.dyna.r.nc`'s `shear_strike`, off by ~2.9e-8 between reruns, well under threshold). Changed the metadata check to compare variable sets and attrs only. Covered by new unit tests in `testsys/unit/test_check_comparisons.py` (within-threshold-but-not-bit-exact -> SUCCESS; differing attrs -> FAIL).
  * Fix - `testsys/e2e/run_e2e.py` dropped the `plotRuptureDynamics` step from the old `testAll.py`'s per-case sequence, so `fault.dyna.r.nc` (one of `check.test.py`'s `fileNameList` comparisons) was never regenerated. Restored: `create.newcase` -> `case.setup` -> `mpirun` -> `plotRuptureDynamics` -> `check.test.py`.
  * Fix - `install-eqdyna.sh` ran `chmod -R 755 scripts` on every build (PROJECT_RULES.md rule 13), flipping the mode bit on all ~60+ tracked files under `scripts/` -- including data files (`*.m`, `*.txt`, `*.mat`) never meant to be executable -- as a side effect of running the e2e test tier. Replaced with `chmod 755` of only the specific entry-point scripts (`case.setup`, `clean.py`, `create.newcase`, `generateFaultInterface`, `plotRuptureDynamics`, `plotSlipAndRPT`), matching their git-tracked exec bits. README's Installation section updated to match (`chmod 755 install-eqdyna.sh`, dropping the recursive `scripts` chmod); resolves pathway_forward.md item 6.
  * Add - `install-eqdyna.sh` to the exec-bit checks in `testsys/regression/test_create_newcase.py` (same guard family as `create.newcase`'s exec bit): it is invoked as `./install-eqdyna.sh` by `testsys/e2e/run_e2e.py`'s fresh-rebuild step, and a lost exec bit there fails a full e2e run the same way.
  * Add - `test.prev/` to `.gitignore` (byproduct of the e2e tier's rule-8 evidence-preservation step; was showing untracked in `git status`, rule 12).
  * Gate: build green; `testsys/run.py unit` 33/33 SUCCESS; `regression` 2/2 SUCCESS; `e2e` 5/5 cases, `check.test.py` 15/15 comparisons SUCCESS. Verified on ubuntu 22.04, gfortran 11.4.0, Open MPI 4.1.1, at commit c4b13e5 + this release's changes.
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

# Environment
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

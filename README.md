# News in 2026
* 20260913 v5.3.6 release notes
  * Fix - PML region-14 cross-axis bug: `src/computePMLDampingVector.f90` and `src/assembleGlobalKU.f90` (formerly `comdampv.f90`) compared the `y` coordinate against `xmax0`/`xmax2` (an x-axis bound) instead of `ymax0`/`ymax2` at all 4 call sites, and left a `z==zmin2` boundary point outside the damped region. Both closed; the shared `pmlRegionDistance` helper in the new `src/func_lib.f90` now carries a single `boundInclusive` convention per caller. Impact was real and physical: `test.drv.a6` peak fault shear traction moved up to ~29 MPa once the correct PML region is damped. `test.reference.results/test.drv.a6/` regenerated from the fixed binary; guarded by new `testsys/regression/test_pml_region_axes.py` (static scan for the exact cross-axis pattern that caused the bug). Resolves pathway_forward.md item 8.
  * Fix - thermal-pressurization `fric_tp_h` (TP shear-zone half-width, the denominator in every thermop kernel) was never assigned in `src/` -- it read as 0.0 in all friclaw==5 runs since dc108b7 (2022-11-28), even though `case.setup` writes the intended 0.02 m (TPV105-3D spec) into `fric(40)`/`on_fault_vars` slot 40, which was never wired to the thermop reads. `src/updateThermalPressurization.f90` (formerly `thermop.f90`) now reads the per-node `fric(40)` value directly. Verified against a clean-room rebuild of v5.2.0 (the SCEC-verified tag): the fixed tree's `test.tpv1053d` output matches v5.2.0 (h=0.02 m); the pre-fix h=0 behavior is categorically different (~1 s rupture time, ~2.9 m slip, vs. the physical result). `test.reference.results/test.tpv1053d/` regenerated from the fixed binary, which doubles as the rule-10 regression test. Resolves pathway_forward.md item 12.
  * Refactor (output-neutral, bit-gated against the two fixes above) - deleted `misc/` (superseded pre-rename originals, -2188 lines, rule 1); converted tabs to 4 spaces across `src/` (EQquasi convention); removed dead commented-out code (57 lines); consolidated 7 duplicated file-existence preambles in `readInputFiles.f90` into `requireInputFile()`; hardened the CI workflow (gate on exit codes, dropped `chmod -R` and `|| true` masking, run `testsys/run.py all`); renamed 6 straggler files/subroutines to verb-phrase names (`comdampv.f90` -> `computePMLDampingVector.f90`, `thermop.f90` -> `updateThermalPressurization.f90`, `mesh4num.f90` -> `countMeshEntities.f90`, and matching subroutine renames), matching the EQquasi naming convention; restructured `globalvar.f90` with a named `FRIC_SLOT_*` fric()-slot map (118 substitutions of literal indices), hoisted `min_norm`/`max_norm` to module init, and deleted 5 orphaned thermal-pressurization globals (including the stale, never-assigned `fric_tp_h`); extracted `src/func_lib.f90` for the shared PML region-distance cascade, `insertFaultInterface`, and `fb`/`vlm` helpers; consolidated `MPI4arn`'s six near-identical sendrecv blocks and `MPI4NodalQuant`'s dispatch ladders (-85 lines); deduplicated Python post-processing (shared `loadFrtData` in `scripts/lib.py`; `plot_on_fault_vars2` merged behind `--with-state`, -180 lines). All output-neutrality verified by clean rebuild + full e2e rerun matching `test.reference.results/` byte-for-byte where the physics did not change.
  * Add - `python/`, a vectorized Python port of EQdyna (NumPy and JAX) covering the slip-weakening, rate-and-state, and thermal-pressurization friction families. Parity against the Fortran solver is documented in `python/README-parity.md`, established against instrumented spike builds (JAX-CPU measured 1.1-2.6x faster than serial Fortran on the cases exercised); the oracle-dump instrumentation and fixtures behind those numbers were deliberately kept spike-only and are not yet reproducible from the committed tree -- tracked as pathway_forward.md item 14, not claimed here as an in-tree regression gate.
  * Add - 37 new unit tests and regression guards (including `testsys/unit/test_lib.py` additions and `testsys/regression/test_pml_region_axes.py`) covering this release's fixes and refactor.
  * Gate: build green; `testsys/run.py unit` SUCCESS; `regression` 3/3 SUCCESS; `e2e` 5/5 cases, `check.test.py` 15/15 comparisons SUCCESS. Verified on ubuntu 22.04, gfortran 11.4.0, Open MPI 4.1.1, netCDF 4.8.1, at commit b7f1647 + this release's changes.
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

# News in 2026
* 20260909 v5.3.4 release notes
  * Fix - uninitialized `thetaPcTmp` in `NewtonRaphson` (src/faulting.f90) for friclaw==5: the copy-back `thetaPc0 = thetaPcTmp` ran unconditionally but `thetaPcTmp` was only assigned inside the friclaw<5 branch, so fric(23)/frt.txt col 22 held a deterministic uninitialized stack value every step. Physics for friclaw==5 is unaffected (theta_pc is not used); frt.txt cols 1-21 are bit-identical before and after. test.tpv1053d golden reference regenerated from the fixed binary and doubles as the regression test.
  * Fix - src/makefile: FC on ubuntu changed from mpif90.mpich to mpif90; mpif90.mpich does not exist on the target system and the ubuntu build was broken.
  * Fix - restored the executable bit on scripts/create.newcase and scripts/generateFaultInterface. Without it, create.newcase was PATH-shadowed by an unrelated project's script of the same name and the test suite could not run from a clean checkout.
  * New - output_src_evol subroutine (src/library_output.f90) writes per-fault-node final slip rate to binary src_evol files, for on-fault state visualization/AI use.
  * Change - output_gm and output_src_evol now run every time step (previously every 10th step via mod(nt,10)==1) in src/driver.f90.
  * Add - PROJECT_RULES.md, a 14-rule project rule book, and pathway_forward.md (formerly docs/PROJECT_STATUS.md), a living status board of open issues, re-check intervals, and a tasks-done log.
  * Add - .gitignore for bin/, __pycache__/, *.mod, *.pyc, scratch/, src/eqdyna, test/.
  * Update - scripts/plotRuptureDynamics: fixed duplicate subplot-axis variable names and added axis labels.
  * Update - scripts/clean.py: also purge src_evol and src* output files.
  * Update - README.md: reworded collaboration and benchmark-performance sections, added DOI links to references.
  * Known issue - src/netcdf_io.f90:76-105 hardcodes fault index 1 instead of the loop variable in on-fault netcdf-input assignment, so multi-fault (ntotft>=2) on-fault netcdf input is not applied correctly. No current test case exercises ntotft>=2. Tracked in pathway_forward.md item 7.
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
chmod -R 755 install-eqdyna.sh scripts
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

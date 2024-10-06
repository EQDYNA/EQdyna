# News in 2024
* 20241006 v5.3.3 release notes:
  * New - verification against new SCEC/USGS Spontaneous Rupture Code Verification benchmarks [TPV36&37](https://strike.scec.org/cvws/tpv36_37docs.html) for 15 deg shallow dipping thrusting. 
  * New - exclusive model setup python script user_defined_params.py for TPV36&37 are under /case_input/test.tpv36 and case_input/test.tpv37, respectively. 
  * Performance: 512 cores are used for 50 m resolution TPV36 on Lonestar6 at TACC using 4 hours and 40 minutes. 
  * New - previous feature of degeneration of hexahedrons (Hughes, 2000) for complex fault geometry is incorporated in the new EQdyna architecture with TPV36&37.
  * New - autotesting workflow is added on GitHub for developers. 
  * New - supporting MacOS (M3 chip tested).
  * Refctor - rename file and subroutine names for clarification.
  * Reference: Hughes, 2000, The Finite Element Method: Linear Static and Dynamic Finite Element Analysis, Dover Publications.  
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

*```EQdyna```* has been extensively verified against benchmark problems from SCEC/USGS Spontaneous Dynamic Rupture Code Verification Excercise.

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
./install-eqdyna.sh -m ubuntu # ubuntu/ls6
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
* test.tpv104
* test.tpv1053d <br/>
[TPV+number is the naming convention of [SCEC/USGS Spontaneous Rupture Code Verification Excercise](https://strike.scec.org/cvws/).] <br/>

For a customized case, please choose the most relevant predefined compset and modify ```user_defined_param.py``` accordingly. <br/>

# Benchmark computational performance and resource
TPV104:   15 seconds simulation time with 0.008 dt and a total of 1875 time steps. It took 40 CPUs to run 24.20 minutes on Lonestar6.  <br/>

# Note
*```EQdyna```* is still under heavy development and comes without any guaranteed functionality. But we hope *```EQdyna```* would be easy to use and we bear this goal in mind when developing it. 

# Developers and collaboration?
We welcome developers and users. If you are interested in developing and collaborations using *```EQdyna```*, please contact Drs. Benchun Duan (bduan@tamu.edu) or Dunyu Liu (dliu@ig.utexas.edu).

# Reference
1. Duan, B. and D.D. Oglesby (2006). Heterogeneous fault stresses from previous earthquakes and the effect on dynamics of parallel strike-slip faults, J. Geophys. Res., 111, B05309, doi:10.1029/2005JB004138.
2. Duan, B. (2012). Dynamic rupture of the 2011 Mw 9.0 Tohoku-Oki earthquake: Roles of a possible subducting seamount, J. Geophys. Res., 117, B05311, doi:10.1029/2011JB009124.
3. Luo, B. and B. Duan (2018). Dynamics of non-planar thrust faults governed by various friction laws, J. Geophys. Res. Solid Earth, 123, https://doi.org/10.1029/2017JB015320.
4. Liu, D. and B. Duan (2018). Scenario Earthquake and Ground‐Motion Simulations in North China: Effects of Heterogeneous Fault Stress and 3D Basin Structure, Bull. Seismol. Soc. Am., 108(4), 2148-2169, doi:10.1785/0120170374.

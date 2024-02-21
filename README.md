# News in 2024
* 20240221 v5.3.2 release notes
  * New - new link for EQdyna. https://github.com/EQDYNA/EQdyna.git
  * New - new organization EQDYNA is created. 
  * Bug - dz for dipping fault. Verified against TPV10.  
  * Bug - avoid int() in Python and FORTRAN, use round() or nint() instead.
  * Add - MATLAB scripts for GM postprocessing. 
  * Add - str1ToFaultAngle and devStrToStrVertRatio for assigning stresses for plastic models.
  * Add - case test.drv.a6.v2 
  * Refactor - move adjustable parameters out of EQdyna. 
  * Refactor - rename file and function names to reflect their intents, for easy search. 
  * Refactor - driver for seismic wave propagation + faulting only. 
  * Refactor - MPI communications for nodal quantities, driver.f90, depreciate PMLwhg.f90 and contm.f90, refactor qdct3.f90, rename qdct3 to ku, offFaultStationSCEC. 
  * Change - Positive dip angles for faults tilting to y+.
  * Change - Use empirical estimate for memory usage. 
  * Update - ubuntu.env.sh

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
* Duan, B. (2010). Role of initial stress rotations in rupture dynamics and ground motion: A case study with implications for the Wenchuan earthquake, J. Geophys. Res. 115.
* Duan, B. (2012). Dynamic rupture of the 2011 Mw 9.0 Tohoku-Oki earthquake: Roles of a possible subducting seamount, J. Geophys. Res. 117
* Duan, B., D. Liu, and A. Yin (2017). Seismic shaking in the North China Basin expected from ruptures of a possible seismic gap, Geophys. Res. Lett. 44 4855-4862.
* Liu, D. and B. Duan (2018). "Scenario Earthquake and Ground‐Motion Simulations in North China: Effects of Heterogeneous Fault Stress and 3D Basin Structure." BSSA.

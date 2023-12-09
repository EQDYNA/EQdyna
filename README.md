# News in 2023
* 20231209 v5.3.1 release notes
  * Python utility generateFaultInterface is created to generate fractal and dipping fault interface;
  * update the test system to compare the numerical residuals between test and reference results; 
  * support a new test test.tpv10, a 60 deg dipping fault;
  * support a new test test.drv.a6, a fractal fault plastic model for ground motion application;
  * refactoring and bug fixes. 

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

# Dependence
*```EQdyna```* requires the following environments and packages: <br/>
  - FORTRAN compiler and MPI
  - Python3
  - netCDF
  
# Installation
A strengh of *```EQdyna```* is that it doesn't need an external mesh generator. <br/> 
Except for some really geometrically complex fault systems, it works well for many applications. <br/> 

*```EQdyna```* now supports Ubuntu and Lonestar6 at TACC. <br/> 

To install it on Ubuntu, 
```
git clone https://github.com/dunyuliu/EQdyna.git
cd EQdyna
chmod 755 install-eqquasi.sh
./install-eqquasi.sh -m ubuntu
```

To install it on Lonestar6 TACC,
```
./install-eqquasi.sh -m ls6
```

To activate environment variable $EQDYNAROOT and add executables to $PATH,
```
source install-eqquasi.sh
```

# Quick Start Guide
After the installation, you just need three steps to run a pre-defined case. <br/>

First, create a new case with the utility *create_newcase*. <br/> 
*create_newcase* takes in two parameters: <br/> 
  (a) case_dir - the case directory where you want to run the job, and <br/>
  (b) compset  - the predefined case (tpv104, tpv1053d, etc) <br/>

[TPV+number is the naming convention of [SCEC/USGS Spontaneous Rupture Code Verification Excercise](https://strike.scec.org/cvws/).] <br/>

```
create_newcase case_dir compset
```
Second, cd into the case_dir, modifiy the *user_defined_params.py* and then run the utility *case.setup*.
```
case.setup
```
Third, run the utility *case.submit*, which will submit the job to the HPC system.
```
case.submit
```

For customized case, please choose the most relevant predefined case and modify the user_defined_param.py accordingly. <br/>

# Test All 
To routinely test the consistency and stability of *```EQdyna```*, <br/>
run the following command:
```
source $EQDYNAROOT/test-all.sh
```
Four pre-defined test cases will be created and run using 4 CPUs, which will take a few minutes.

# Currently supported compset (more to come)
* drv.a6, for determinisitc ground motion with fractal fault and plasticity
* tpv8
* tpv10
* tpv104
* tpv1053d
* meng2023a, for stress inversion
* meng2023cb

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

# Past release notes
* 20231130 v5.3.0 release notes
  * refactor faulting.f90 and meshgen.f90;
  * introduce a system for quick testing by running testAll.py; 
  * test system now supports tpv8, tpv104, tpv1053d, tpv1053d.6c, meng2023a, meng2023cb;
  * update input parameter system with a single defaultParameters.py and customized user_defined_params.py;

* 20230327 
  * *```EQdyna```* works on Ubuntu now. 
  * install-eqdyna.sh can support multiple systems. 
  * test-all.sh will create and run four pre-defined cases with 4 CPUs with a few minutes. 
* 20230321 *```eqdyna.docker```* is published via [dunyuliu/eqdyna.docker](https://hub.docker.com/repository/docker/dunyuliu/eqdyna.docker/general) 

# News in 2022
* A set of python utlities that make *```EQdyna```* much easier to use.
* Thermal pressurization implemented and benchmarked against TPV105-3D.
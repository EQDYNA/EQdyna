# News in 2022
* A set of python utlities that make *```EQdyna```* much easier to use.
* Thermal pressurization implemented and benchmarked against TPV105-3D.

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
*```EQdyna```* requires FORTRAN and MPI compiler, python, and netCDF libaries. Currently, some post-processing scripts are written in MATLAB but those will be phased out to python scripts.

# Installation on Lonestar6 at TACC, UT Austin.
A strengh of *```EQdyna```* is that it doesn't need an external mesh generator. Except for some really geometrically complex fault systems, it generally works for lots of applications. So, its installation is pretty easy and straightforward. Although install-eqdyna.sh is written for Lonestar6, what *```EQdyna```* needs are generally available from HPCs. If the netCDF libaries is missing, it could cost you some time to install and link it to *```EQdyna```*.   
```
git clone https://github.com/dunyuliu/EQdyna.git
cd EQdyna
source install-eqdyna.sh
```

# Quick Start Guide
After the installation, you just need three steps to run a predefined case. <br/>

First, create a new case with the utility *create_newcase*. <br/> 
*create_newcase* takes in two parameters: <br/> 
  (a) case_dir - the case directory where you want to run the job, and <br/>
  (b) compset  - the predefined case (tpv104, tpv105, etc) <br/>

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

# Example
A good starting example would be compset==tpv104 (benchmark problem 104, dynamic rupture with rate- and state- friction with strong rate weakening friction law, 100 m resolution, 40 CPUs, a few minutes for 3 seconds in simulation time.)

# Note
*```EQdyna```* is still under heavy development and comes without any guaranteed functionality. But we hope *```EQdyna```* would be easy to use and we bear this goal in mind when developing it. 

# Developers and collaboration?
We welcome developers and users. If you are interested in developing and collaborations using *```EQdyna```*, please contact Drs. Benchun Duan (bduan@tamu.edu) or Dunyu Liu (dliu@ig.utexas.edu).

# Reference
* Duan, B. (2010). Role of initial stress rotations in rupture dynamics and ground motion: A case study with implications for the Wenchuan earthquake, J. Geophys. Res. 115.
* Duan, B. (2012). Dynamic rupture of the 2011 Mw 9.0 Tohoku-Oki earthquake: Roles of a possible subducting seamount, J. Geophys. Res. 117
* Duan, B., D. Liu, and A. Yin (2017). Seismic shaking in the North China Basin expected from ruptures of a possible seismic gap, Geophys. Res. Lett. 44 4855-4862.
* Liu, D. and B. Duan (2018). "Scenario Earthquake and Ground‐Motion Simulations in North China: Effects of Heterogeneous Fault Stress and 3D Basin Structure." BSSA.

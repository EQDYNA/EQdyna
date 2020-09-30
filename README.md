# EQdyna 3D

EQdyna 3D is a parallel finite element software to simulate earthquake spontaneous dynamic 
ruptures and seismic wave propagations for geometrically complex fault systems written in 
FORTRAN90. EQdyna 3D is highly efficient and scalable due to its adoptions of explicit 
time integration, underintegrated hexahedral elements, degenerated wedge elements, hourglass 
controls, perfectly matched layer and parallelization. Frequency independent Q by coarsed
-grained memory scheme allows seismic attenuation in a time marching FEM. Frictional 
constitutations such as the classical slip-weakening and various forms of rate- and state-
friction enables applications to dynamic ruptures, ground motion, fully dynamic earthquake
cycles, etc. Drucker-prager off-fault viscoplasticity allows integration of numerical models 
with near-fault geologic observations. EQdyna has been verified against benchmark problems
from community-led SCEC/USGS code verification excercise. Recently, termalpressurization is 
implemented.

### Author:  Dunyu Liu and Bin Luo
### Date:    09/30/2020
### Contact: dunyuliu@tamu.edu
## Version 5.1.0; Git tag v5.1.0; Parent 5.0.0.
# Major changes:
* Thermal pressurization (tp) is implemented;
* Library.f90 is created to store non-major subroutines;
* Input file FE_Fric.txt is added and accordingly subroutine readfric in Read_Input_Files.f90;
* Additional dimensions in array fric for tp;
* Input file FE_Stations.txt is added and on- and off-fault stations can be flexibly assinged;
* [Verification] This version of code has been verified in the benchmark problem SCEC TPV105 3D;
* Detailed and other minor changes please refer to changelog.md.

### Dunyu Liu, 08/12/2020
# Version 5.0.0
## Features
* A major structural change occurs: controllable parameters, model and fault information, 
	material properties are all moved to input .txt files. This signals the developing goals 
	to move adjustable quantities out of the bones of finite element calculations.
* Read_Input_Files.f90 is added to load input 4 files FE_Global.txt, FE_Model_Geometry.txt, 
	FE_Model_Geometry.txt, FE_Material.txt. 
* The code runs smoothly and is verified against TPV104 in terms of rupture time contour.
* FE_Fric.txt is planned as an input but not developed in this version.

### Dunyu Liu, Bin Luo, 10/04/2016, dunyuliu@gmail.com
# Version 4.2
## Features
* Add the rate- and state- friction with aging law (friclaw=3) and slip law (friclaw=4) 
	(Bin Luo) incorporated by Dunyu.
* Time expense analysis with timeanalysis.m.
* This version is verified against SCEC TPV104.

### Dunyu Liu, 10/04/2016, dunyuliu@gmail.com
# Version 4.1.2
## Patch update
* Test a plastic model.
* Disable switches formma and formkd in driver.f90.
* Add several parameters related to plastic model in globalvar.f90

### Dunyu Liu, 10/04/2016, dunyuliu@gmail.com
# Version 4.1.1
## Features
* Significant simplification of the system.
* Reduction in variable declarations.
* Moving all controllable parameters into globalvar.f90.
* This version is tested against SCEC TPV8.

### Dunyu Liu, Bin Luo, 09/29/2016, dunyuliu@gmail.com 
# Version 4.0
## Features
* 3D MPI (Bin Luo) incorporated with Version 3.2.1 by Dunyu Liu.
* Further simplification.
* More controllable parameters in globalvar.f90.
* This version is verified against SCEC TPV 8.
* EOS retired and batch file written for ADA. 

### Dunyu Liu, 09/19/2015, dunyuliu@gmail.com
# Version 3.2.1
## Features
* Coarse-grained Q model is implemented (Ma and Liu, 2006; Day, 1998).
* Elastic and plastic models are combined. 
* qconstant.f90 is added.
* Controllable parameters are moved to globalvar.f90 and switched to change mechanicms are added.
* The code is validated against the model with PML and Q model in Ma and Liu (2006) and is later used in the Tianjin Scenario Earthquake project (Duan et al., 2017; Liu and Duan, 2018).

### Dunyu Liu, 09/2015, dunyuliu@gmail.com
# Version 3.1.2
## Patch update
* Double-couple point source is implemented based on Ma and Liu (2006)

### Dunyu Liu, 09/2015, dunyuliu@gmail.com
# Version 3.1.1
## Features
* Implementation of the Perfectly Matched Layer (Ma and Liu, 2006)
* Major update to use 1D array to store nodal forces, and kinematic quantities. It is designed to solve the different dimensions of regular and PML elements. It involves many files that pass such quantities. 
* Element type is introduced to sperate regular and PML elements.
* 2 F90 files are added and 13 are modified. 
* Both viscous and KF78 hourglass controls are implemented.
* Makefile updated to include the 2 new files. 
* This version is verified against SCEC TPV29&30. 

### Benchun Duan, 09/12/2014, bduan@tamu.edu
# Version 3.1
## Features
* Finite Element Method (FEM) is based on Duan (2010)
* Traction-at-split Node (TSN) technique in faulting.f90 to simulate earthquake faults (Day et al., 2005)
* Drucker-Prager plastic yielding (Ma and Andrews, 2010).
* MPI/OpenMP hybrid parallelization (Wu et al., 2011)
* The code is verified in SCEC benchmark problems 
  TPV18-21 (scecdata.usc.edu/cvws/) and used in Jingqian Kang's PhD project - low velocity fault zone response to nearby earthquake. 

## Model description
* Material properties, initial stresses, friction, and plasticity (c = 0 Mpa case) of Ma and Andrews (2010) are used.
* Revised to taper stress drop at two lateral edges of the fault by following the scheme for tapering stresses at top and bottom parts as they proposed, to gradually stop rupture.
* The fault dimension is x=[-16.1 km, 16.1 km] and 
  z=[-15.1 km, 0], with element size of 100 m.

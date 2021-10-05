# EQdyna

EQdyna is a parallel Finite Element software to simulate earthquake spontaneous dynamic rupture and seismic wave propagation. It aims to simulate earthquakes on geometrically complex fault systems. It is written in FORTRAN90.
EQdyna is highly efficient and scalable due to its adoptions of explicit time integration, underintegrated hexahedral elements, degenerated wedge elements, hourglass controls, perfectly matched layer and parallelization with MPI. Frequency independent Q by coarsed-grained memory scheme allows seismic attenuation in a time marching FEM. Frictional 
constitutions such as slip-weakening and various forms of rate- and state-friction enables applications to dynamic ruptures, ground motion, fully dynamic earthquake cycles, etc. Drucker-prager off-fault viscoplasticity allows integration of numerical models with near-fault geologic observations. EQdyna has been verified against benchmark problems from community-led SCEC/USGS code verification excercise. Recently, termal pressurization is 
implemented. EQdyna is part of the fully dynamic earthquake cycle simulator EQsimu (Liu et al., 2020, GJI).

# Software requirements
* a Unix-lie operating system
* build tools cmake
* Fortran compilers
* a function MPI environment
* Matlab installed to use the post-processing script

# Note
EQdyna is still under heavy development and comes without any guaranteed functionality. At the moment, we cannot provide much support for general users.

# Collaboration
If you are interested in collaboration using EQdyna, please contact Drs. Benchun Duan (bduan@tamu.edu) or Dunyu Liu (dliu@ig.utexas.edu).

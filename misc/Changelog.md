# Version 5.2.0
* new file:   LICENSE
* modified:   README.md
* modified:   code/PMLwhg.f90
* modified:   code/Read_Input_Files.f90
* modified:   code/comdampv.f90
* modified:   code/contm.f90
* modified:   code/driver.f90
* modified:   code/eqdyna3d.f90
* modified:   code/faulting.f90
* deleted:    code/formlm.f90
* modified:   code/fric.f90
* modified:   code/globalvar.f90
* modified:   code/hrglss.f90
* new file:   code/library_degeneration.f90
* new file:   code/library_output.f90
* modified:   code/makefile
* modified:   code/mesh4num.f90
* modified:   code/meshgen.f90
* deleted:    code/parcon.f90
* modified:   code/qdckd.f90
* modified:   code/qdcshg.f90
* modified:   code/qdcshl.f90
* deleted:    code/qdct1.f90
* modified:   code/qdct2.f90
* modified:   code/qdct3.f90
* modified:   code/thermop.f90
* deleted:    code/thermop_new.f90
* deleted:    code/thermop_new2.f90
* deleted:    code/thermop_old.f90
* new file:   code/warning.f90
* renamed:    script/rtp3DMPI.m -> script/Plot_rpt.m
* deleted:    Guide_EQdyna_4.2.md
* deleted:    make/makefile
# Detailed changes:
New: a license is attached. 
New: warning.f90 is introduced. It alerts users about potential limitations 
	in EQdyna.
New: library_output.f90 is created to handle output. Subroutines include 
	output_onfault_st, output_offfault_st, output_frt, output_timeanalysis. 
	Future output could be implemented here.
New: library_degeneration.f90 is created. 
	It contains currently five subroutines: wedge, wedge4num, tetra, tetra4num,
	and reorder. Wedge and tetra are used in meshgen after creating an element 
	to decide whether to do degeneration. 
	
	A domain is determined in wedge and tetra to be one element off the fault. 
	[Note: The region could be specified or be related to realistic fault geometry, in future development.] 
	wedge4num and tetra4num are used in mesh4num to create additional number 
	counts in numel due to degeneration. 
	
	The current function is to degenerate a hexahedron to two wedges or five tetrahedrons. 
	
	Subroutine reorder is created to simplify the process of degeneration by 
	assigning the order of local node number in an element. Repeated indices 
	indicate collapsing of nodes for degeneration. 
	[Note: the shape of wedge/tetrahedrons are important. The edges of equivalent edges are desired.]  
New: a controlling parameter C_degen is introduced. 
	0 stands for no degeneration. 
	1 stands for degeneration to wedges.
	2 stands for degeneration into tetrahedrons. 
New: mesh quality of wedges are improved in eqdyna when C_degen == 1.
Improvement: move arrays that don’t change over simulations, i.e., global arrays, 
	to globalvar.f90. No array transference in major subroutines.
Improvement: adopt real kind = dp, dp = selected_kind(15,307). 
Improvement: initials are assigned to double precision.
Improvement: hourglass control is applied to every element regardless of element types. 
Improvement: Shape functions are kept constant for all types of elements. 
	No need to sum up local shape functions for collapsed nodes in degeneration.
Verification: C_degen = 0, 1, 2 for TPV105 3D in SCEC/USGS code verification project.   

# Version 5.1.0
* modified:   Changelog.md
* modified:   README.md
* renamed:    batch/batchADA3DMPI.txt -> batch/ADArun.txt
* new file:   batch/TERRArun.slurm
* deleted:    batch/run256eos.hyb
* deleted:    batch/runvanerne.inta
* modified:   code/Read_Input_Files.f90
* modified:   code/driver.f90
* modified:   code/eqdyna3d.f90
* modified:   code/faulting.f90
* modified:   code/globalvar.f90
* new file:   code/library.f90
* new file:   code/makefile
* modified:   code/meshgen.f90
* new file:   code/thermop.f90
* new file:   code/thermop_new.f90
* new file:   code/thermop_new2.f90
* new file:   code/thermop_old.f90
* deleted:    code/vlm.f90
* new file:   document/Guide_EQdyna_4.2.md
* modified:   input/FE_Fault_Geometry.txt
* new file:   input/FE_Fric.txt
* modified:   input/FE_Global.txt
* modified:   input/FE_Material.txt
* modified:   input/FE_Model_Geometry.txt
* new file:   input/FE_Stations.txt
* new file:   script/rm.tcsh
* new file:   script/rtp3DMPI.m
* deleted:    Guide_EQdyna_4.2.md
* deleted:    make/makefile
# Detailed changes:
New feature: thermal pressurization is done with Eq (17) in Andrews (2002), 
	which calculates the pore pressure change. Eq (14) can be readily modified 
	to calculate temperature change due to tp. Use dt in the elastodynamic 
	process to do the temporal integration. The estimation is good even though
	a much smaller timestep would be needed to accurately calculate temperature
	and pore pressure evolution with a finite difference method. 
New feature: Library.f90 is created to store miscellaneous functions and 
	subroutines. Subroutine vlm is moved to library.f90. 
New feature: functions B1, B2 and B3 in TPV105-3D to taper parameters are added 
	in Library.f90  
New feature: FE_Fric.txt is created to store inputs of frictional parameters 
	and FE_Fric.txt is read in by subroutine readfric. A new friction law type
	friclaw == 5 is introduced to account for thermal pressurization. 
New feature: subroutine readfric is added in Read_Input_Files.f90.
New feature: Subroutine thermop is created. It is to calculate temperature and
	pore pressure evolution on a newly created mesh, which is in a much finer 
	scale compared to the FE mesh and close to the fault. The pore pressure and 
	temperature at integer time step are stored in array fric’s 51th and 52th 
	dimensions, respectively. The pore pressure is transferred to subroutine
	faulting when updating the effective normal stress tnrm. 
[not used]: numtp and dxtp are introduced and assigned in subroutine globalvar. 
	Numtp stands for the number of the dimension of tp mesh that is perpendicular
	to the fault surface. Dxtp is the spatial resolution of the tp mesh in that
	direction. They are used to allocate arrays in subroutine thermop and to 
	calculate second derivatives of T and p relative to y coordinates. 
New feature: define a set of new parameters with fric_ in their names and they 
	are declared as global variables in subroutine globalvar. They are read in 
	by calling subroutine readfric.
Improvement: the array fric is expanded to 100 in its first dimensions. It is 
	expanded more than needed for later development. 
New feature: a new input file FE_Stations.txt is introduced. It is read in by 
	subroutines readstations1 (for station numbers for both on-fault and on-
	surface ones) and readstations2 (for their coordinates).  
New feature: subroutines readstations1 and readstations2 are added in Read_Input_Files.f90. 
Verification: The results are verified in the benchmark problem TPV105-3D in 
	SCEC/USGS code verification project.   

# Version 5.0.0
* modified:   README.md
* new file:   code/Read_Input_Files.f90
* modified:   code/eqdyna3d.f90
* modified:   code/globalvar.f90
* modified:   code/meshgen.f90
* new file:   code/runTERRA.slurm
* new file:   input/FE_Fault_Geometry.txt
* new file:   input/FE_Global.txt
* new file:   input/FE_Material.txt
* new file:   input/FE_Model_Geometry.txt
* modified:   make/makefile

# Version 4.2

## Modified
* driver.f90
* eqdyna3d.f90
* faulting.f90
* globalvar.f90
* meshgen.f90
* fric.f90

# Version 4.1.2

## Modified
* driver.f90
* eqdyna3d.f90
* faulting.f90
* globalvar.f90
* meshgen.f90
* qdckd.f90
* qdct3.f90

# Version 4.1.1

## Modified
* PMLwhg.f90
* comdampv.f90
* driver.f90
* eqdyna3d.f90
* faulting.f90
* globalvar.f90
* hrglss.f90
* mesh4num.f90
* meshgen.f90
* qdckd.f90
* qdct2.f90
* qdct3.f90
* makefile

# Version 4.0

## Modified
* PMLwhg.f90
* comdampv.f90
* driver.f90
* eqdyna3d.f90
* faulting.f90
* globalvar.f90
* hrglss.f90
* mesh4num.f90
* meshgen.f90
* parcon.f90
* qdckd.f90
* qdct2.f90
* qdct3.f90
* makefile

## Added 
* batchADA3DMPI.txt

# Version 3.2.1

## Modified
* PMLwhg.f90
* comdampv.f90
* driver.f90
* eqdyna3d.f90
* faulting.f90
* fric.f90
* globalvar.f90
* hrglss.f90
* mesh4num.f90
* meshgen.f90
* parcon.f90
* qdckd.f90
* qdct3.f90
* makefile

## Added
* qconstant.f90

# Version 3.1.2

## Modified 
* driver.f90
* faulting.f90
* mesh4num.f90
* meshgen.f90
* parcon.f90
* qdckd.f90

# Version 3.1.1 

## Modified: 
* contm.f90
* driver.f90
* eqdyna3d.f90
* faulting.f90
* formlm.f90
* hrglss.f90
* mesh4num.f90
* meshgen.f90
* parcon.f90
* prop3d.f90
* qdckd.f90
* qdct2.f90
* qdct3.f90

* makefile

## Added
* PMLwhg.f90
* comdampv.f90


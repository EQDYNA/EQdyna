#!/bin/bash

# Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
# MIT

# Environment variable MACHINE is defined in install-eqdyna.sh.
# Currently, EQdyna supports: 
#     ls6
#     ubuntu
#     grace
SYSTEM = ${MACHINE}

ifeq ($(SYSTEM), ubuntu)
    FFLAGS     = -fopenmp -ffree-line-length-none -O3
    FC         = mpif90
    LIB        = /usr/lib/x86_64-linux-gnu
    INC        = /usr/include
    NETCDF_LIB = -L${LIB} -lnetcdf -lnetcdff
    NETCDF_INC = -I${INC}
    
else ifeq ($(SYSTEM), ls6)
    FFLAGS     = -qopenmp
    FC         = mpiifort
    LIB        = ${TACC_MKL_LIB}
    INC        = ${TACC_MUMPS_INC}
    NETCDF_INC = -I${TACC_NETCDF_INC}
    NETCDF_LIB = -L${TACC_NETCDF_LIB} -lnetcdf -L${TACC_NETCDF_LIB} -lnetcdff 
    
else ifeq ($(SYSTEM), grace)
    FFLAGS     = 
    FC         = mpiifort
    LIB        = ${EBROOTNETCDF}/lib64
    INC        = ${EBROOTNETCDF}/include
    NETCDF_INC = -I${INC}
    NETCDF_LIB = -L${LIB} -lnetcdf -L${LIB} -lnetcdff 

else ifeq ($(SYSTEM), macos)
    FFLAGS     = -O3 -ffast-math -ftree-vectorize 
    FC         = mpif90
    LIB        = ${MACOS_NETCDF_LIB}
    INC        = ${MACOS_NETCDF_INC}
    NETCDF_INC = -I${MACOS_NETCDF_INC}
    NETCDF_LIB = -L${MACOS_NETCDF_LIB} -lnetcdf -L${MACOS_NETCDFF_LIB} -lnetcdff 
endif

OPT  = -c ${FFLAGS} 
OPT1 = ${FFLAGS}

OBJ  = eqdyna3d.o  globalvar.o countMeshEntities.o  meshgen.o   calcLocalShapeFunc.o    \
       driver.o    assembleGlobalMass.o     calcGlobalShapeFunc.o                calcB.o      \
       assembleGlobalKU.o        calcElemMass.o    calcElemKU.o     calcHourglassResist.o    fric.o      \
       faulting.o              computePMLDampingVector.o  calcQAttenuationCoeff.o             \
       readInputFiles.o        updateThermalPressurization.o   library.o   checkInputConsistency.o   \
       library_degeneration.f90            library_output.o        \
       netcdf_io.o

clean:eqdyna
	rm -r *.o

eqdyna: $(OBJ) 
	$(FC) $(OPT1) $(OBJ) -o eqdyna ${NETCDF_LIB}
eqdyna3d.o: eqdyna3d.f90 globalvar.o driver.o countMeshEntities.o \
            meshgen.o readInputFiles.o checkInputConsistency.o         \
            library_output.o netcdf_io.o
	$(FC) $(OPT) eqdyna3d.f90
globalvar.o: globalvar.f90
	$(FC) $(OPT) globalvar.f90
countMeshEntities.o: countMeshEntities.f90 globalvar.o library_degeneration.o
	$(FC) $(OPT) countMeshEntities.f90
meshgen.o: meshgen.f90 globalvar.o library.o library_degeneration.o
	$(FC) $(OPT) meshgen.f90
calcLocalShapeFunc.o: calcLocalShapeFunc.f90 globalvar.o
	$(FC) $(OPT) calcLocalShapeFunc.f90		
driver.o: driver.f90 assembleGlobalMass.o assembleGlobalKU.o calcHourglassResist.o faulting.o \
          globalvar.o computePMLDampingVector.o updateThermalPressurization.o
	$(FC) $(OPT) driver.f90
assembleGlobalMass.o: assembleGlobalMass.f90 calcGlobalShapeFunc.o library.o globalvar.o
	$(FC) $(OPT) assembleGlobalMass.f90
calcGlobalShapeFunc.o: calcGlobalShapeFunc.f90 globalvar.o
	$(FC) $(OPT) calcGlobalShapeFunc.f90
library.o: library.f90 globalvar.o
	$(FC) $(OPT) library.f90
assembleGlobalKU.o: assembleGlobalKU.f90 calcElemKU.o globalvar.o 
	$(FC) $(OPT) assembleGlobalKU.f90
calcElemMass.o: calcElemMass.f90 globalvar.o
	$(FC) $(OPT) calcElemMass.f90
calcElemKU.o: calcElemKU.f90 calcB.o calcQAttenuationCoeff.o globalvar.o
	$(FC) $(OPT) calcElemKU.f90
calcB.o: calcB.f90 globalvar.o
	$(FC) $(OPT) calcB.f90
calcHourglassResist.o: calcHourglassResist.f90 globalvar.o
	$(FC) $(OPT) calcHourglassResist.f90
fric.o: fric.f90 globalvar.o
	$(FC) $(OPT) fric.f90
faulting.o: faulting.f90 fric.o globalvar.o
	$(FC) $(OPT) faulting.f90
computePMLDampingVector.o: computePMLDampingVector.f90
	$(FC) $(OPT) computePMLDampingVector.f90
calcQAttenuationCoeff.o: calcQAttenuationCoeff.f90
	$(FC) $(OPT) calcQAttenuationCoeff.f90
readInputFiles.o: readInputFiles.f90 globalvar.o
	$(FC) $(OPT) readInputFiles.f90
updateThermalPressurization.o: updateThermalPressurization.f90 globalvar.o
	$(FC) $(OPT) updateThermalPressurization.f90
checkInputConsistency.o: checkInputConsistency.f90 globalvar.o
	$(FC) $(OPT) checkInputConsistency.f90
library_degeneration.o: library_degeneration.f90 globalvar.o 
	$(FC) $(OPT) library_degeneration.f90
library_output.o: library_output.f90 globalvar.o 
	$(FC) $(OPT) library_output.f90
netcdf_io.o: netcdf_io.f90 globalvar.o 
	$(FC) $(OPT) netcdf_io.f90 ${NETCDF_INC}

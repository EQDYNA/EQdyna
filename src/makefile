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

# PYDUMP=1 builds a distinct, non-default `eqdyna-pydump` binary that
# generates the Python-port parity fixtures (testsys/parity/). driver.f90
# and eqdyna3d.f90 call pydump_step/pydump_state UNCONDITIONALLY (no
# preprocessor, no #ifdef) -- the seam is at LINK time, not compile time:
# the default build links pydump_noop.o (empty subroutines, call+return),
# PYDUMP=1 links pydump.o (the real dump writer). This keeps eqdyna3d.o and
# driver.o's compiled bytes IDENTICAL between the two builds (verified:
# `cmp` on both objects after independent compiles). An earlier `#ifdef
# PYDUMP` + `-cpp` design was tried and rejected: enabling `-cpp`, even
# with the macro inactive, measurably changed the compiled object's SIZE
# (not just embedded metadata) for reasons not fully isolated in the time
# available -- so it was dropped in favor of this cleaner link-time seam,
# which does not require `-cpp` anywhere. Residual, non-eliminable
# difference from the true pre-pydump pristine source: two extra
# unconditional CALL instructions in eqdyna3d.o/driver.o (into an empty
# subroutine in the default build) -- see testsys/parity/README.md for the
# functional-equivalence verification (e2e case comparison) that stands in
# for byte-identity against pristine.
PYDUMP ?= 0
ifeq ($(PYDUMP), 1)
    PYDUMP_OBJ = pydump.o
else
    PYDUMP_OBJ = pydump_noop.o
endif

OPT  = -c ${FFLAGS}
OPT1 = ${FFLAGS}

OBJ  = eqdyna3d.o  globalvar.o errorCodes.o countMeshEntities.o  meshgen.o   calcLocalShapeFunc.o    \
       driver.o    assembleGlobalMass.o     calcGlobalShapeFunc.o                calcB.o      \
       assembleGlobalKU.o        calcElemMass.o    calcElemKU.o     calcHourglassResist.o    fric.o      \
       faulting.o              computePMLDampingVector.o  calcQAttenuationCoeff.o             \
       readInputFiles.o        updateThermalPressurization.o   library.o   func_lib.o   checkInputConsistency.o   \
       library_degeneration.o              library_output.o        \
       netcdf_io.o             $(PYDUMP_OBJ)

# `clean:eqdyna` (no space) made clean DEPEND on eqdyna, so `make clean` built
# the binary before deleting the objects. Nothing in the repo invokes it, which
# is why it went unnoticed. Also removes .mod files, which were never cleaned.
.PHONY: clean
clean:
	rm -f *.o *.mod

eqdyna: $(OBJ)
	$(FC) $(OPT1) $(OBJ) -o eqdyna ${NETCDF_LIB}
eqdyna-pydump:
	$(MAKE) PYDUMP=1 eqdyna
	mv eqdyna eqdyna-pydump
eqdyna3d.o: eqdyna3d.f90 globalvar.o driver.o countMeshEntities.o errorCodes.o \
            meshgen.o readInputFiles.o checkInputConsistency.o         \
            library_output.o netcdf_io.o $(PYDUMP_OBJ)
	$(FC) $(OPT) eqdyna3d.f90
pydump.o: pydump.f90 globalvar.o
	$(FC) $(OPT) pydump.f90
pydump_noop.o: pydump_noop.f90
	$(FC) $(OPT) pydump_noop.f90
globalvar.o: globalvar.f90
	$(FC) $(OPT) globalvar.f90
errorCodes.o: errorCodes.f90
	$(FC) $(OPT) errorCodes.f90
countMeshEntities.o: countMeshEntities.f90 globalvar.o library_degeneration.o
	$(FC) $(OPT) countMeshEntities.f90
meshgen.o: meshgen.f90 globalvar.o func_lib.o library_degeneration.o errorCodes.o
	$(FC) $(OPT) meshgen.f90
calcLocalShapeFunc.o: calcLocalShapeFunc.f90 globalvar.o
	$(FC) $(OPT) calcLocalShapeFunc.f90		
driver.o: driver.f90 assembleGlobalMass.o assembleGlobalKU.o calcHourglassResist.o faulting.o errorCodes.o \
          globalvar.o computePMLDampingVector.o updateThermalPressurization.o
	$(FC) $(OPT) driver.f90
assembleGlobalMass.o: assembleGlobalMass.f90 calcGlobalShapeFunc.o library.o func_lib.o globalvar.o errorCodes.o
	$(FC) $(OPT) assembleGlobalMass.f90
calcGlobalShapeFunc.o: calcGlobalShapeFunc.f90 globalvar.o errorCodes.o
	$(FC) $(OPT) calcGlobalShapeFunc.f90
library.o: library.f90 globalvar.o
	$(FC) $(OPT) library.f90
func_lib.o: func_lib.f90 globalvar.o errorCodes.o
	$(FC) $(OPT) func_lib.f90
assembleGlobalKU.o: assembleGlobalKU.f90 calcElemKU.o globalvar.o func_lib.o
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
computePMLDampingVector.o: computePMLDampingVector.f90 globalvar.o func_lib.o errorCodes.o
	$(FC) $(OPT) computePMLDampingVector.f90
calcQAttenuationCoeff.o: calcQAttenuationCoeff.f90
	$(FC) $(OPT) calcQAttenuationCoeff.f90
readInputFiles.o: readInputFiles.f90 globalvar.o errorCodes.o
	$(FC) $(OPT) readInputFiles.f90
updateThermalPressurization.o: updateThermalPressurization.f90 globalvar.o
	$(FC) $(OPT) updateThermalPressurization.f90
checkInputConsistency.o: checkInputConsistency.f90 globalvar.o errorCodes.o
	$(FC) $(OPT) checkInputConsistency.f90
library_degeneration.o: library_degeneration.f90 globalvar.o 
	$(FC) $(OPT) library_degeneration.f90
library_output.o: library_output.f90 globalvar.o 
	$(FC) $(OPT) library_output.f90
netcdf_io.o: netcdf_io.f90 globalvar.o errorCodes.o
	$(FC) $(OPT) netcdf_io.f90 ${NETCDF_INC}

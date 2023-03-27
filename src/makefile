#!/bin/bash

# This makefile file is part of software EQdyna.

# This makefile is to compile EQdyna v5.2.3. 

# Environment variable MACHINE is defined in install-eqdyna.sh.
# Currently, EQdyna supports: 
#     ls6
#     ubuntu
SYSTEM = ${MACHINE}

ifeq ($(SYSTEM), ubuntu)
    FFLAGS     = -fopenmp -ffree-line-length-none
    FC         = mpif90.mpich
    LIB        = /usr/lib/x86_64-linux-gnu
    INC        = /usr/include
    NETCDF_LIB = -L${LIB} -lnetcdf -L{LIB} -lnetcdff
    NETCDF_INC = -I${INC}
    
else ifeq ($(SYSTEM), ls6)
    FFLAGS     = -qopenmp
    FC         = mpiifort
    LIB        = ${TACC_MKL_LIB}
    INC        = ${TACC_MUMPS_INC}
    NETCDF_INC = -I${TACC_NETCDF_INC}
    NETCDF_LIB = -L${TACC_NETCDF_LIB} -lnetcdf -L${TACC_NETCDF_LIB} -lnetcdff 
    
endif

OPT  = -c ${FFLAGS} 
OPT1 = ${FFLAGS}

OBJ  = eqdyna3d.o  globalvar.o mesh4num.o  meshgen.o   qdcshl.o    \
       driver.o    qdct2.o     qdcshg.o    contm.o     qdcb.o      \
       qdct3.o     contma.o    qdckd.o     hrglss.o    fric.o      \
       faulting.o  PMLwhg.o    comdampv.o  qconstant.o             \
       readInputFiles.o        thermop.o   library.o   warning.o   \
       library_degeneration.f90            library_output.o        \
       netcdf_io.o

clean:eqdyna
	rm -r *.o

eqdyna: $(OBJ) 
	$(FC) $(OPT1) $(OBJ) -o eqdyna ${NETCDF_LIB}
eqdyna3d.o: eqdyna3d.f90 globalvar.o driver.o mesh4num.o \
            meshgen.o readInputFiles.o warning.o         \
            library_output.o netcdf_io.o
	$(FC) $(OPT) eqdyna3d.f90 
globalvar.o: globalvar.f90
	$(FC) $(OPT) globalvar.f90 	
mesh4num.o: mesh4num.f90 globalvar.o library_degeneration.o
	$(FC) $(OPT) mesh4num.f90
meshgen.o: meshgen.f90 globalvar.o library.o library_degeneration.o
	$(FC) $(OPT) meshgen.f90
qdcshl.o: qdcshl.f90 globalvar.o
	$(FC) $(OPT) qdcshl.f90		
driver.o: driver.f90 qdct2.o qdct3.o hrglss.o faulting.o \
          globalvar.o comdampv.o thermop.o
	$(FC) $(OPT) driver.f90
qdct2.o: qdct2.f90 qdcshg.o contm.o library.o globalvar.o
	$(FC) $(OPT) qdct2.f90
qdcshg.o: qdcshg.f90 globalvar.o
	$(FC) $(OPT) qdcshg.f90
contm.o: contm.f90 globalvar.o
	$(FC) $(OPT) contm.f90
library.o: library.f90 globalvar.o
	$(FC) $(OPT) library.f90
qdct3.o: qdct3.f90 contma.o qdckd.o globalvar.o PMLwhg.o
	$(FC) $(OPT) qdct3.f90
contma.o: contma.f90 globalvar.o
	$(FC) $(OPT) contma.f90
qdckd.o: qdckd.f90 qdcb.o qconstant.o globalvar.o
	$(FC) $(OPT) qdckd.f90
qdcb.o: qdcb.f90 globalvar.o
	$(FC) $(OPT) qdcb.f90
hrglss.o: hrglss.f90 globalvar.o
	$(FC) $(OPT) hrglss.f90
fric.o: fric.f90 globalvar.o
	$(FC) $(OPT) fric.f90
faulting.o: faulting.f90 fric.o globalvar.o
	$(FC) $(OPT) faulting.f90
PMLwhg.o: PMLwhg.f90 globalvar.o
	$(FC) $(OPT) PMLwhg.f90 
comdampv.o: comdampv.f90
	$(FC) $(OPT) comdampv.f90
qconstant.o: qconstant.f90
	$(FC) $(OPT) qconstant.f90
readInputFiles.o: readInputFiles.f90 globalvar.o
	$(FC) $(OPT) readInputFiles.f90 
thermop.o: thermop.f90 globalvar.o
	$(FC) $(OPT) thermop.f90 
warning.o: warning.f90 globalvar.o
	$(FC) $(OPT) warning.f90
library_degeneration.o: library_degeneration.f90 globalvar.o 
	$(FC) $(OPT) library_degeneration.f90
library_output.o: library_output.f90 globalvar.o 
	$(FC) $(OPT) library_output.f90
netcdf_io.o: netcdf_io.f90 globalvar.o 
	$(FC) $(OPT) netcdf_io.f90 ${NETCDF_INC}

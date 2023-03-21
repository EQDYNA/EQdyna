# On LS6 at TACC, UTA.
netcdf_incdir = ${TACC_NETCDF_INC}
netcdf_libdir = ${TACC_NETCDF_LIB}

F90   = ${FC}
OPT1  = 
OPT   = -c -ffree-line-length-none

OBJ   = eqdyna3d.o globalvar.o mesh4num.o meshgen.o qdcshl.o driver.o qdct2.o qdcshg.o contm.o qdcb.o \
		qdct3.o contma.o qdckd.o hrglss.o fric.o faulting.o PMLwhg.o comdampv.o qconstant.o \
		readInputFiles.o thermop.o library.o warning.o library_degeneration.f90 \
		library_output.o netcdf_io.o
		
clean:eqdyna
	rm -r *.o
eqdyna: $(OBJ) 
	$(F90) $(OPT1) $(OBJ) -o eqdyna -L$(netcdf_libdir) -lnetcdf -L$(netcdf_libdir) -lnetcdff
eqdyna3d.o: eqdyna3d.f90 globalvar.o driver.o mesh4num.o meshgen.o readInputFiles.o warning.o library_output.o netcdf_io.o
	$(F90) $(OPT) eqdyna3d.f90 
globalvar.o: globalvar.f90
	$(F90) $(OPT) globalvar.f90 	
mesh4num.o: mesh4num.f90 globalvar.o library_degeneration.o
	$(F90) $(OPT) mesh4num.f90
meshgen.o: meshgen.f90 globalvar.o library.o library_degeneration.o
	$(F90) $(OPT) meshgen.f90
qdcshl.o: qdcshl.f90 globalvar.o
	$(F90) $(OPT) qdcshl.f90		
driver.o: driver.f90 qdct2.o qdct3.o hrglss.o faulting.o globalvar.o comdampv.o thermop.o
	$(F90) $(OPT) driver.f90
qdct2.o: qdct2.f90 qdcshg.o contm.o library.o globalvar.o
	$(F90) $(OPT) qdct2.f90
qdcshg.o: qdcshg.f90 globalvar.o
	$(F90) $(OPT) qdcshg.f90
contm.o: contm.f90 globalvar.o
	$(F90) $(OPT) contm.f90
library.o: library.f90 globalvar.o
	$(F90) $(OPT) library.f90		
qdct3.o: qdct3.f90 contma.o qdckd.o globalvar.o PMLwhg.o
	$(F90) $(OPT) qdct3.f90
contma.o: contma.f90 globalvar.o
	$(F90) $(OPT) contma.f90
qdckd.o: qdckd.f90 qdcb.o qconstant.o globalvar.o
	$(F90) $(OPT) qdckd.f90
qdcb.o: qdcb.f90 globalvar.o
	$(F90) $(OPT) qdcb.f90
hrglss.o: hrglss.f90 globalvar.o
	$(F90) $(OPT) hrglss.f90	
fric.o: fric.f90 globalvar.o
	$(F90) $(OPT) fric.f90	
faulting.o: faulting.f90 fric.o globalvar.o
	$(F90) $(OPT) faulting.f90	
PMLwhg.o: PMLwhg.f90 globalvar.o
	$(F90) $(OPT) PMLwhg.f90 
comdampv.o: comdampv.f90
	$(F90) $(OPT) comdampv.f90
qconstant.o: qconstant.f90
	$(F90) $(OPT) qconstant.f90
readInputFiles.o: readInputFiles.f90 globalvar.o
	$(F90) $(OPT) readInputFiles.f90 
thermop.o: thermop.f90 globalvar.o
	$(F90) $(OPT) thermop.f90 
warning.o: warning.f90 globalvar.o
	$(F90) $(OPT) warning.f90
library_degeneration.o: library_degeneration.f90 globalvar.o 
	$(F90) $(OPT) library_degeneration.f90
library_output.o: library_output.f90 globalvar.o 
	$(F90) $(OPT) library_output.f90
netcdf_io.o: netcdf_io.f90 globalvar.o 
	$(F90) $(OPT) netcdf_io.f90 -I${netcdf_incdir}	

F90 = mpiifort
OPT2=-ipo -O3 -xHost -openmp -heap-arrays -limf -c 
OPT=-c
OPT1= 
OBJ = eqdyna3d.o globalvar.o mesh4num.o \
         meshgen.o qdcshl.o driver.o \
	 qdct2.o qdcshg.o contm.o qdcb.o vlm.o qdct3.o \
	 contma.o qdckd.o hrglss.o fric.o faulting.o \
	 PMLwhg.o comdampv.o qconstant.o Read_Input_Files.o
eqdyna3drun: $(OBJ) 
	$(F90) $(OPT1) $(OBJ) -o eqdyna-hyb
eqdyna3d.o: eqdyna3d.f90 globalvar.o driver.o mesh4num.o meshgen.o Read_Input_Files.o
	$(F90) $(OPT) eqdyna3d.f90 
globalvar.o: globalvar.f90
	$(F90) $(OPT) globalvar.f90 	
mesh4num.o: mesh4num.f90 globalvar.o
	$(F90) $(OPT) mesh4num.f90
meshgen.o: meshgen.f90 globalvar.o
	$(F90) $(OPT) meshgen.f90
qdcshl.o: qdcshl.f90 globalvar.o
	$(F90) $(OPT) qdcshl.f90		
driver.o: driver.f90 qdct2.o qdct3.o hrglss.o faulting.o globalvar.o comdampv.o
	$(F90) $(OPT) driver.f90
qdct2.o: qdct2.f90 qdcshg.o contm.o vlm.o globalvar.o
	$(F90) $(OPT) qdct2.f90
qdcshg.o: qdcshg.f90 globalvar.o
	$(F90) $(OPT) qdcshg.f90
contm.o: contm.f90 globalvar.o
	$(F90) $(OPT) contm.f90
vlm.o: vlm.f90 globalvar.o
	$(F90) $(OPT) vlm.f90		
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
Read_Input_Files.o: Read_Input_Files.f90 globalvar.o
	$(F90) $(OPT) Read_Input_Files.f90 
clean:
	rm fslipout_*
	rm fsurout_*
	rm frt*
	rm fnode*

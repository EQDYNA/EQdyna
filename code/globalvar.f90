!...To MPI parallel, comment definitions for large arrays.
! B.D. 12/20/08
!
MODULE globalvar
  implicit none
  save  
  !...inout/output units
  integer (kind=4) :: iingl=11,ioutgl=12, ioutft=13,ioutrt1=14,&
  	   ioutrt2=15,ioutst=16,ioutsl=17,ioutgm=18,ioutoff=19,&
           ioutrat=20,ioutdp=21
  !...measure portion CPU time
  real (kind=8), dimension (7) :: timeused=0.0
  real (kind=8) :: time1, time2, time=0.0, btime=0.0,term
  external gethrtime	!this function in the SUN Fortran library
  !...execution control parameters	   	   
  integer (kind=4) :: iexec,irank,numseq,nsd,ndof,numeg
  !...scalar variables (eps. for hexahedral elements)
  integer (kind=4),parameter :: nen=8,ned=3,nee=24,nesd=3,nrowsh=4, &
           nrowb=6,nrowc=6,nstr=6,noid=2
  real (kind=8) :: w	!integration weight	   	   
  !...definitions for faults
  integer (kind=4),parameter :: ntotft=1
  integer (kind=4) :: nucfault,friclaw
  real (kind=8) :: critd0, &
  		mus,mud,cohes,brangle,tv
  !...time sequence parameters
  real (kind=8) :: dt	!actual time step 
  integer (kind=4) :: nplpts,nstep,nhplt,nhplt1,nhshw
  !...time histories	   	    
  integer (kind=4) :: locplt,myrec=0
!-------------------------------------------------------------------!
!-----------------------Sep.19.2015/ D.Liu--------------------------!
!-------------Controllable parameters for EQdyna V3.2.1-------------! 
integer(kind=4)::C_elastic=1
!	1=elastic version;
!	0=plastic version;
integer(kind=4)::C_Q=0!
!AAA:Only works with C_elastic==1.
!	1=allow Q attenuation;
!	0=do not allow Q.
integer(kind=4)::C_hg=1
!	1=KF78 hourglass control(HG);
!	2=Viscous HG
integer(kind=4)::C_dc=0!In driver.f90
!AAA: fault should be fixed.(mus==10000.)
!	1=allow double couple(DC) point source;
!	0=do not allow.
integer(kind=4)::nPML=6
!	Thickness (counted by nodes) of PML.
real(kind=8)::R=0.01! Theoretical reflection coefficient for PML
!	Other options for pairs of (nPML/R)
!	nPML/R=6/0.01;10/0.001;20/0.0001.Collino& Tsogka(2001)
real(kind=8)::kapa_hg=0.1!Coefficient for viscous HG.
!AAA:Only works when C_hg==2.
!	Typical valus varies from 0.05~0.15.Goudreau& Hallquist(1982).
real(kind=8) :: rat=1.0
!	Enlarge ratio for buffers to use 
integer(kind=4)::numat=2!Material types.
real(kind=8)::dx!Element size. In parcon.f90.
real(kind=8)::critt0=0.4,srcrad0=4000.,vrupt0=1500.
real(kind=8)::xsource=-12.2e3,ysource=0.0,zsource=-13.8e3
integer(kind=4)::ninterval=1
real(kind=8)::surxmax=0e3,surxmin=0e3,surymax=0e3,surymin=0e3
real(kind=8)::vmaxPML=5716.
!	For 3D MPI Partitioning
integer(kind=4)::npx=6,npy=6,npz=2
logical,dimension(6)::fltMPI
integer(kind=4),dimension(6)::fltnum=0
integer(kind=4),allocatable,dimension(:)::fltgm
integer(kind=4),allocatable,dimension(:)::fltl,fltr,fltf,fltb,fltd,fltu
integer (kind=4),dimension(9)::numcount
!--------End of Controllable parameters for EQdyna V3.2.1-----------!  
!-------------------------------------------------------------------!    
end MODULE globalvar 	      

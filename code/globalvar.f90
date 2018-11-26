!...To MPI parallel, comment definitions for large arrays.
! B.D. 12/20/08
!
MODULE globalvar
  implicit none
  save  
  !
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
           nrowb=6,nrowc=6,nstr=6,noid=2,numat=1
  real (kind=8) :: w	!integration weight	   	   
  !...definitions for faults
  integer (kind=4),parameter :: ntotft=1
  integer (kind=4) :: nucfault,friclaw
  real (kind=8) :: critd0,critt0,srcrad0,vrupt0,xsource,ysource,zsource, &
  		mus,mud,cohes,brangle,tv
  !...time sequence parameters
  real (kind=8) :: dt	!actual time step 
  integer (kind=4) :: nplpts,nstep,nhplt,nhplt1,nhshw
  !...time histories	   	    
  integer (kind=4) :: locplt,myrec=0
  !...global arrays to dea lwith possible zero arrays. B.D. 8/12/10
!  integer (kind=4),allocatable,dimension(:,:) :: idhist
!  real (kind=8),allocatable,dimension(:,:) :: dout
  !...global fault arrays to deal with zero fault node for some processor. 
  !only allocate once in meshgen.f90 whrn fault node is nonzero. B.D. 10/17/09
!  integer (kind=4),allocatable,dimension(:,:) :: anonfs
!  integer (kind=4),allocatable,dimension(:,:,:) :: nsmp
!  real (kind=8),allocatable,dimension(:,:) :: fnft,arn,r4nuc,arn4m,slp4fri
!  real (kind=8),allocatable,dimension(:,:,:) :: fric,un,us,ud
!  real (kind=8),allocatable,dimension(:,:,:) :: fltsta
	     
end MODULE globalvar 	      

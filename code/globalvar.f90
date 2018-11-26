!===================================================================!
MODULE globalvar
implicit none
save  
integer(kind=4)::iingl=11,ioutgl=12, ioutft=13,ioutrt1=14,&
  	ioutrt2=15,ioutst=16,ioutsl=17,ioutgm=18,ioutoff=19,&
    ioutrat=20,ioutdp=21
real(kind=8),dimension(7)::timeused=0.0
real(kind=8)::time1,time2,time=0.0, btime=0.0
external gethrtime	!this function in the SUN Fortran library   	   
integer(kind=4)::iexec=1,irank=0,numseq=1,nsd=3,ndof=3,numeg=1
integer(kind=4)::nen=8,ned=3,nee=24,nesd=3,nrowsh=4,nrowb=6,nrowc=6,nstr=6,noid=2
real(kind=8)::w	!integration weight	   	   
real(kind=8)::critd0,cohes,brangle
integer(kind=4)::nplpts,nstep,nhplt=1,nhplt1=2,nhshw=1   	    
integer(kind=4)::locplt=1,myrec=0
real(kind=8)::pi=3.1415926535897931,rdampm=0.0,rdampk=0.1
!-------------------------------------------------------------------!
!-------------Controllable parameters for EQdyna3d V4.0-------------! 
integer(kind=4)::C_elastic=0
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
real(kind=8)::rat=1.025
!	Enlarge ratio for buffers to use 
real(kind=8)::dx=100.,dy,dz!Element size.
integer(kind=4)::dis4uniF=200,dis4uniB=200
!Number of uniform sized elements normal to the fault
real(kind=8)::rhow=1000.,b11=0.926793,b33=1.073206,b13=-0.169029,&
	critt0=0.5,srcrad0=4000.,vrupt0=1500.,&!Info on forced rupture.
	bulk=0.1934,coheplas=1.36e6,tv=0.03,ccosphi,sinphi,mus,mud
real(kind=8)::xsource=-5e3,ysource=0.0,zsource=-10e3!Nucleation point.
!===================================================================!
!Specify informations on on- and off- fault stations and
! model and fault geometries. 
integer(kind=4)::ninterval=1,nftmx,nonmx,nonfs(1)=12,n4nds=18,an4nds(2,18)
real(kind=8)::xonfs(2,12,1),x4nds(3,18)
real(kind=8)::surxmax=0e3,surxmin=0e3,surymax=0e3,surymin=0e3
real(kind=8)::xmin=-40e3,xmax=40e3,ymin=-40e3,ymax=35e3,zmin=-50e3,zmax=0.0
real(kind=8)::fxmin=-20e3,fxmax=20e3,fymin=0.0,fymax=0.0,fzmin=-20e3,fzmax=0.0
real(kind=8)::fstrike=270.,fdip=90.
real(kind=8)::fltxyz(2,4,1)
character(len=15)::projectname='SCECTPV27v2',author='D.Liu'
!===================================================================!
!Specify maximum Vp for PML and timing information
real(kind=8)::vmaxPML=6000.,term=20.,dt=0.008
integer(kind=4)::nucfault=1,friclaw=1,ntotft=1
real(kind=8),allocatable,dimension(:,:)::dout
integer(kind=4),allocatable,dimension(:,:)::idhist
!===================================================================!
!3DMPI Partitioning along x/y/z axis, repect.
integer(kind=4)::npx=5,npy=5,npz=4,master=0,nprocs
logical,dimension(6)::fltMPI
integer(kind=4),dimension(6)::fltnum=0
integer(kind=4),allocatable,dimension(:)::fltgm
integer(kind=4),allocatable,dimension(:)::fltl,fltr,fltf,fltb,fltd,fltu
integer (kind=4),dimension(9)::numcount
!-------------------------------------------------------------------!    
end MODULE globalvar 	      

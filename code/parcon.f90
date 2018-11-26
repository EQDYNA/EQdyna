SUBROUTINE parcon(xmin,xmax,ymin,ymax,zmin,zmax,fltxyz,&
  emod,pois,mushr,rho,vp,vs,rdampm,rdampk,ccosphi,sinphi)
!...provide control and property parameters for 3D FEM mesh and dynamic 
!	simulation. B.D. 4/21/09
!...incorporate into EQdyna. B.D. 10/28/09
!...add in material properties calculations, e.g., calculating elastic 
!  constants from given Vp, Vs, and rho. B.D. 8/20/10
!...add in Drucker-Prager plasticity parameter. B.D. 1/6/12
!
  use globalvar
  implicit none
  real (kind=8)::pi=3.1415926535897931
  !...define non-global parameters
  integer(kind=4)::i,j
  real(kind=8)::fltwdt,fltx1,fltx2,flty,fltz1,fltz2, &
  	xmin,xmax,ymin,ymax,zmin,zmax,temp
  real (kind=8),dimension(numat) :: emod,pois,mushr,rho,rdampm,rdampk, &
       vp,vs,lamda,coh4p,tanphi,ccosphi,sinphi
  real (kind=8),dimension(2,4,ntotft) :: fltxyz !4:x,y,z,strike/dip
  !  	
  !...control parameters that seldom need to change.
  iexec=1
  irank=0
  numseq=1
  nsd=3
  ndof=3
  numeg=1
  !
  !...termination, output, & simulation time info
  term=15.0
  dt=0.008
  nhplt = 1
  nhplt1 = 2
  nhshw = 1
  nstep = idnint(term/dt)
  !
  !...material properties: normally Vp, Vs are given, we need to 
  !  calculate emod and pois here. B.D. 8/20/10
  !  similarly, given cohesion, and internal friction, calculate
  !  c*cos(phi) and sin(phi). B.D. 1/6/12
!For Double-couple point source. Ma and Liu (2006)   
  ! rho=(/2600.,2700./)
  ! vp=(/2800.,6000./)
  ! vs=(/1500.,3464./)
!For TPV8  
rho=(/2700.,2700./)
  vp=(/5716.,5716./)
  vs=(/3300.,3300./)
  coh4p=(/0.e6,0.0/)
  tanphi=(/0.75,0.75/)
  do i=1,numat
    mushr(i) = rho(i) * vs(i) * vs(i)
    lamda(i) = vp(i) * vp(i) * rho(i) - 2. * mushr(i)
    emod(i) = mushr(i)*(3*lamda(i)+2*mushr(i))/(lamda(i)+mushr(i))
    pois(i) = 0.5 * lamda(i) / (lamda(i)+mushr(i))
    temp = atan(tanphi(i))
    ccosphi(i) = coh4p(i) * dcos(temp)
    sinphi(i) = dsin(temp)
  enddo
  rdampm=(/0.0,0.0/)
  rdampk=(/0.1,0.1/)
  do i=1,numat
    rdampk(i) = rdampk(i) * dt	!related to dt
  enddo
  !
  !...fault and friction
  nucfault=1
  friclaw=1
  !critd0=0.5
  tv= 0.03 !time for viscoplasticity: s wave travel one element size
  critt0=0.08  !time-weakening for nucleation
  srcrad0=4000.0
  vrupt0=2000.0
  xsource=0.0
  ysource=0.0
  zsource=-5000.
  !mus=0.76
  !mud=0.448
  !cohes=0.2e6
  !
!For Double-couple point source. Ma and Liu (2006)    
  ! fltxyz=reshape((/0,100,0,0,-100,0,270,90/),&
    ! (/2,4,1/))
!For TPV8	
  fltxyz=reshape((/-15e3,15e3,0,0,-15e3,0,270,90/),&
    (/2,4,1/))
  !...convert degree to arc for strike and dip of faults
  do i=1,ntotft
    do j=1,2
      fltxyz(j,4,i)=fltxyz(j,4,i) * pi / 180.
    enddo
  enddo
  !...branch angle in arc.
!  brangle = fltxyz(1,4,2) - fltxyz(1,4,1)
  !strike=270.
  !dip=90.
  !
  !fltwdt=19575.
  !fltx1=-24075.
  !fltx2=24075.
  !flty=0.0
  !fltz1=-19575.
  !fltz2=0.0
  !
  !...element size and model boundary
!For Double-couple point source. Ma and Liu (2006)  
	! dx=100
	! xmin=-10e3
	! xmax=10e3
	! ymin=-10e3
	! ymax=10e3
	! zmin=-12e3
	! zmax=0
	dx=100
	xmin=-40e3
	xmax=40e3
	ymin=-40e3
	ymax=40e3
	zmin=-50e3
	zmax=0	
end SUBROUTINE parcon

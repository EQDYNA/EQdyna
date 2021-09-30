!/* Copyright (C) 2006-2020, Earthquake Modeling Lab @ Texas A&M University. 
! * All Rights Reserved.
! * This code is part of software EQdyna, please see EQdyna License Agreement
! * attached before you copy, download, install or use EQdyna./
MODULE globalvar

	implicit none
	save

	character(len=30) :: mm, sttmp, dptmp, bodytmp, projectname='SCECTPV105-3D', &
						author='D.Liu & B.Luo'	
	character(len=90) :: loca
	integer, parameter :: dp = selected_real_kind(15,307)
	real (kind = dp) :: timeused(9)=0.0d0, time1, time2, time=0.0d0, btime=0.0d0, &
						timebegin, timeover 	
	integer (kind = 8) :: time_array(8) 	
	!-------------------------------------------------------------------!
	!------------------           Constants              ---------------! 
	real (kind = dp) :: pi=3.1415926535897931d0, rdampm=0.0d0, rdampk=0.1d0, grav = 9.8d0
	integer (kind = 4) :: np = 1000000
	!-------------------------------------------------------------------!
	!------------------            FE system             ---------------!
	integer (kind=4) :: nsd=3, ndof=3, nen=8, ned=3, nee=24, nesd=3, nrowsh=4, &
						nrowb=6, nrowc=6, nstr=6, noid=2, numnp, numel, neq, maxm, maxs, &
						npx, npy, npz, master=0, me, nprocs						
	real (kind = dp) :: w, term, dt, totmemcost, memcost=0.0d0
	real (kind = dp), allocatable, dimension(:) :: brhs, alhs, v1, d1, s1, miuonf, vponf, eleporep, pstrain, eledet, fnms 
	real (kind = dp), allocatable, dimension(:,:) :: x, d, v, mat, shl, fnft, arn, r4nuc, arn4m, slp4fri, state, elemass, ss, plane1, plane2
	real (kind = dp), allocatable, dimension(:,:,:) :: fric, un, us, ud, fltsta, fltslp, eleshp, phi
	integer (kind = 4), allocatable, dimension(:) :: nftnd,  id1, ids, et, locid, dof1, surface_node_id
	integer (kind = 4), allocatable, dimension(:,:) :: ien, anonfs
	integer (kind = 4), allocatable, dimension(:,:,:) :: nsmp
	!-------------------------------------------------------------------!
	!------------------     Geometry and Material        ---------------!
	real (kind = dp) :: rat, dx, dy, dz, xsource, ysource, zsource, &
						xmin, xmax, ymin, ymax, zmin, zmax, &
						xmin1, xmax1, ymin1, ymax1, zmin1, zmax1
	real (kind = dp), allocatable, dimension(:) :: fxmin, fxmax, fymin, fymax, fzmin, fzmax
	real (kind = dp), allocatable, dimension(:,:) :: material
	real (kind = dp), allocatable, dimension(:,:,:) :: fltxyz						
	integer (kind = 4) :: dis4uniF, dis4uniB, nmat, n2mat	
	
	integer(kind=4) :: nplpts, nstep, nhplt=1, locplt=1
	!-------------------------------------------------------------------!
	!------------------               PML                ---------------!
	integer (kind = 4) :: nPML = 6 !Thickness (counted by nodes) of PML.nPML=6
	real (kind = dp) :: R=0.01d0! Theoretical reflection coefficient for PML
	!	Other options for pairs of (nPML/R)
	!	nPML/R=6/0.01;10/0.001;20/0.0001.Collino& Tsogka(2001)
	real (kind = dp) :: PMLb(8), vmaxPML=6000.0d0
	!-------------------------------------------------------------------!
	!-------------         Controllable parameters         -------------! 
	integer(kind=4)::C_elastic! 1=elastic version; 0=plastic version;
	integer(kind=4)::C_nuclea !1=allow artificial nucleation; !0=disabled;
	integer(kind=4)::C_degen ! 0=brick; 1=wedge; 2=textra
	integer(kind=4)::nucfault,friclaw,ntotft
	integer(kind=4)::C_Q=0! !AAA:Only works with C_elastic==1.
	! 1=allow Q attenuation; 0=do not allow Q.
	integer(kind=4)::C_hg=1 !1=KF78 hourglass control(HG); 2=Viscous HG
	integer(kind=4)::C_dc=0!In driver.f90
	!AAA: fault should be fixed.(mus==10000.)
	!	1=allow double couple(DC) point source;
	!	0=do not allow.
	real(kind = dp)::kapa_hg=0.1d0!Coefficient for viscous HG.
	!AAA:Only works when C_hg==2.
	!	Typical valus varies from 0.05~0.15.Goudreau& Hallquist(1982).
	real(kind = dp)::rhow=1000.d0, b11=0.926793, b33=1.073206, b13=-0.169029, &
		critt0=0.5d0, srcrad0=4000.d0, vrupt0=1500.d0, critd0, cohes, brangle, &
		bulk=0.75d0,coheplas=5.0d6, tv, ccosphi, sinphi, mus, mud
	integer(kind=4)::ninterval=1,nftmx,nonmx, nt
	real(kind = dp)::fstrike=270.,fdip=90.
	
	logical,dimension(6)::fltMPI
	integer(kind=4),dimension(6)::fltnum=0
	integer(kind=4),allocatable,dimension(:)::fltgm
	integer(kind=4),allocatable,dimension(:)::fltl,fltr,fltf,fltb,fltd,fltu
	integer (kind=4),dimension(9)::numcount
	!===================================================================!
	!Thermop, TPV105 3D
	real(kind = dp), allocatable, dimension(:,:) :: Tatnode, patnode
	real(kind = dp) :: dxtp,tpw! m
	!Parameters for friction laws
	real(kind = dp):: fric_sw_fs, fric_sw_fd, fric_sw_D0, &
	fric_rsf_a, fric_rsf_deltaa0, fric_rsf_b, fric_rsf_Dc, fric_rsf_r0, fric_rsf_v0, &
	fric_rsf_vinix, fric_rsf_viniz, fric_rsf_fw, fric_rsf_vw, fric_rsf_deltavw0, &
	fric_tp_a_th, fric_tp_rouc, fric_tp_lambda, fric_tp_h, &
	fric_tp_a_hy, fric_tp_deltaa_hy0, fric_ww, fric_w, fric_ini_sliprate,&
	fric_tp_pini, fric_tp_Tini
	real(kind = dp), allocatable, dimension(:,:,:,:):: frichis
	!-------------------------------------------------------------------!
	!-------------                   Output                -------------! 
	real(kind = dp),allocatable,dimension(:,:)::dout
	integer(kind=4),allocatable,dimension(:,:)::idhist	
	! integer (kind = 4) :: n4nds = 6, an4nds(2,6), nonfs(1) = 13 ! on-surface station nubmer
	! real (kind = dp) :: xonfs(2,13,1), x4nds(3,6) 
	integer (kind = 4) :: n4nds, n4onf, n4out, ndout ! on-surface station nubmer
	integer (kind = 4), allocatable :: an4nds(:,:), nonfs(:), n4yn(:)
	!an4nds(2,n4nds): for indexing on-surface stations in different MPI prcocs.
	!nonfs(ntotft): numbers of on-fault stations for multiple faults.
	real (kind = dp), allocatable :: xonfs(:,:,:), x4nds(:,:)
	!xonfs(2,nonfs(ntotft),ntotft)
	!x4nfs(3,n4nds)
	logical :: lstr
	integer(kind = 4) :: TPV = -1, output_plastic
	!For arrays of rough fault geometry
	real (kind = dp) :: dxtmp, perturb = 0.7d0, gamar = 0.66d0, roumax = 2.8d3
	integer (kind = 4) :: nnx, nnz, rough_fault = 1, timeinfo = 0, output_ground_motion = 1, surface_nnode = 0
	real (kind = dp), allocatable, dimension(:,:) :: rough_geo	
	!-------------------------------------------------------------------!    
end MODULE globalvar 	      

!/* Copyright (C) 2006-2023, Earthquake Modeling Lab @ Texas A&M University. 
! * All Rights Reserved.
! * This code is part of software EQdyna, please see EQdyna License Agreement
! * attached before you copy, download, install or use EQdyna./
MODULE globalvar

    implicit none
    save

    character(len=30) :: mm,    &
        sttmp,  dptmp,  bodytmp,    projectname,    author
    character(len=90) :: loca
    
    logical :: lstr, fltMPI(6)
    
    integer, parameter :: dp = selected_real_kind(15,307)
    
    real (kind = dp) :: timeused(9) = 0.0d0,&
        time1,      time2,      time=0.0d0, &
        btime=0.0d0,            timebegin,  timeover,   pi=4*atan(1.0_dp), &
        rdampm=0.0d0,           rdampk=0.1d0,           grav = 9.8d0, &
        w,          term,       dt,         totmemcost, memcost=0.0d0,&
        rat,        dx,         dy,         dz,         xsource,    &
        ysource,    zsource,    xmin,       xmax,       ymin,       &
        ymax,       zmin,       zmax,       xmin1,      xmax1,      &
        ymin1,      ymax1,      zmin1,      zmax1,      nucR,       &
        nucRuptVel, nucdtau0,   max_norm,   min_norm,   &
        ! R: Theoretical reflection coefficient for PML.
        ! Other options for pairs of (nPML/R) are:
        !   nPML/R=6/0.01;10/0.001;20/0.0001.Collino& Tsogka(2001)
        R = 0.01d0, &
        PMLb(8),    vmaxPML=6000.0d0, &
        ! kapa_hg: Coefficient for viscous HG.
        ! Typical valus varies from 0.05~0.15. Goudreau& Hallquist(1982).
        kapa_hg = 0.1d0,    &
        rhow=1000.d0,           b11=0.926793,           b33=1.073206, &
        b13=-0.169029,          critt0=0.2d0,           srcrad0=2500.d0,&
        vrupt0=1500.d0,         critd0,     cohes,      brangle, &
        bulk=0.75d0,            coheplas=5.0d6,         tv,     &
        ccosphi,    sinphi,     mus,        mud,        fstrike=270.d0, &
        fdip=90.d0, dxtp,       tpw, &
        fric_sw_fs, fric_sw_fd, fric_sw_D0, fric_rsf_a, fric_rsf_deltaa0,&
        fric_rsf_b, fric_rsf_Dc,            fric_rsf_r0,                &
        fric_rsf_v0,            fric_rsf_vinix,         fric_rsf_viniz, &
        fric_rsf_fw,            fric_rsf_vw,            fric_rsf_deltavw0, &
        fric_tp_a_th,           fric_tp_rouc,           fric_tp_lambda, &
        fric_tp_h,  fric_tp_a_hy,           fric_tp_deltaa_hy0,         &
        fric_ww, fric_w,        fric_ini_sliprate,      fric_tp_pini,   &
        fric_tp_Tini,           dxtmp,      perturb = 0.7d0,            &
        gamar = 0.66d0,         roumax = 2.8d3,         rough_fx_min,   &
        rough_fx_max,           rough_fz_min
        
    integer (kind = 8) :: time_array(8)

    integer (kind = 4)  :: np = 1000000, &
        nsd=3,  ndof=3, nen=8,  ned=3,  nee=24, nesd=3, nrowsh=4, &
        nrowb=6,nrowc=6,nstr=6, noid=2, numnp,  numel,  neq,      & 
        maxm,   maxs,   npx,    npy,    npz,    master=0,         &
        me,     nprocs, mode, &
        nplpts, nstep,  nhplt=1,        locplt=1, &
        dis4uniF,       dis4uniB,       nmat,   n2mat,  ninterval=1,&
        nftmx,  nonmx,  nt,     TPV = -1,       output_plastic,     &
        nnx,    nnz,    rough_fault,    timeinfo = 0,   outputGroundMotion,&
        surface_nnode = 0,  &
        ! Thickness (counted by nodes) of PML.nPML=6
        nPML = 6, &
        n4nds,  n4onf,  n4out,  ndout, &
        ! 1=elastic version; 0=plastic version;
        C_elastic,  &
        ! 1=allow artificial nucleation; !0=disabled;
        C_nuclea,   &
        ! 0=brick; 1=wedge; 2=textra
        C_degen,    &
        nucfault,   friclaw,    ntotft, &
        ! Only works with C_elastic==1.
        ! 1=allow Q attenuation; 0=do not allow Q.
        C_Q = 0,    &
        !C_hg: 1, KF78 hourglass control(HG); 2, Viscous HG
        C_hg = 1,   &
        ! C_dc: double coutple source; 1, yes; 0, no.
        C_dc = 0,   &
        fltnum(6) = 0,  &
        numcount(9)
     
    real (kind = dp), allocatable, dimension(:) :: brhs,   alhs,&
        v1,     d1,     s1,     miuonf, vponf,  eleporep,       &
        pstrain,eledet, fnms,   fxmin,  fxmax,  fymin,  fymax,  &
        fzmin,  fzmax
    real (kind = dp), allocatable, dimension(:,:) :: x, &
        d,      v,      mat,    shl,    fnft,   arn,    r4nuc,  &
        arn4m,  slp4fri,state,  elemass,ss,     plane1, plane2, &
        Tatnode,patnode,dout,   material,       rough_geo,      &
        x4nds
    real (kind = dp), allocatable, dimension(:,:,:) :: fric,    &
        un,     us,     ud,     fltsta, fltslp, eleshp, phi,    &
        fltxyz, xonfs
    real (kind = dp), allocatable, dimension(:,:,:,:):: frichis

    integer (kind = 4), allocatable, dimension(:) :: nftnd,     &
        id1,    ids,    et,     locid,  dof1,   surface_node_id,&
        nonfs,  n4yn,   fltl,   fltr,   fltf,   fltb,   fltd,   &
        fltu,   fltgm
    integer (kind = 4), allocatable, dimension(:,:) :: ien,     &
        anonfs, idhist, an4nds
    integer (kind = 4), allocatable, dimension(:,:,:) :: nsmp

end MODULE globalvar           

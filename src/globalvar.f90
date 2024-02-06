! Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
! MIT
MODULE globalvar

    implicit none
    character(len=30) :: mm,    &
        sttmp,  dptmp,  bodytmp,    projectname,    author
    character(len=90) :: stLocStamp
    logical :: fltMPI(6)
    integer, parameter :: dp=selected_real_kind(15,307)
    real (kind = dp) :: compTimeInSeconds(9)=0.0d0, &
        startTimeStamp,         timeElapsed=0.0d0, &
        MPICommTimeInSeconds=0.0d0,            simuStartTime,  pi=4*atan(1.0_dp), &
        rdampm=0.0d0,           rdampk,                 grav = 9.8d0, &
        w,          totalSimuTime,       dt,         totmemcost, memcost=0.0d0,&
        rat,        dx,         dy,         dz,         xsource,    &
        ysource,    zsource,    xmin,       xmax,       ymin,       &
        ymax,       zmin,       zmax,       xmin1,      xmax1,      &
        ymin1,      ymax1,      zmin1,      zmax1,      nucR,       &
        nucRuptVel, nucdtau0,   max_norm,   min_norm,   &
        ! R: Theoretical reflection coefficient for PML.
        ! Other options for pairs of (nPML/R) are:
        !   nPML/R=6/0.01;10/0.001;20/0.0001.Collino& Tsogka(2001)
        R = 0.01d0, &
        PMLb(8),    vmaxPML, &
        ! kapa_hg: Coefficient for viscous HG.
        ! Typical valus varies from 0.05~0.15. Goudreau& Hallquist(1982).
        kapa_hg = 0.1d0,    &
        tol     = 1.0d-5,   &
        rhow,                   critt0=0.2d0,           srcrad0=2500.d0,&
        vrupt0=1500.d0,         critd0,     cohes,      brangle, &
        bulk,                   coheplas,               tv,     &
        ccosphi,    sinphi,     fstrike,        &
        fdip,       slipRateThres, nucT, &
        fric_sw_fs, fric_sw_fd, fric_sw_D0, fric_rsf_a, fric_rsf_deltaa0,&
        fric_rsf_b, fric_rsf_Dc,            fric_rsf_r0,                &
        fric_rsf_v0,            fric_rsf_vinix,         fric_rsf_viniz, &
        fric_rsf_fw,            fric_rsf_vw,            fric_rsf_deltavw0, &
        fric_tp_a_th,           fric_tp_rouc,           fric_tp_lambda, &
        fric_tp_h,  fric_tp_a_hy,           fric_tp_deltaa_hy0,         &
        fric_ww, fric_w,        fric_ini_sliprate,      fric_tp_pini,   &
        fric_tp_Tini,           dxtmp,            &
        gamar,                  roumax,         rough_fx_min,   &
        rough_fx_max,           rough_fz_min,  &
        str1ToFaultAngle,       devStrToStrVertRatio
        
    integer (kind = 8) :: dateTimeStamp(8)

    integer (kind = 4)  :: np = 1000000, &
        nsd=3,  ndof=3, nen=8,  ned=3,  nee=24, nesd=3, nrowsh=4, &
        nrowb=6,nrowc=6,nstr=6, totalNumOfNodes,  totalNumOfElements,  totalNumOfEquations,      & 
        sizeOfEqNumIndexArr,   sizeOfStressDofIndexArr,   npx,    npy,    npz,    masterProcsId=0, &
        me,     totalNumOfMPIProcs, mode,   nstep, &
        dis4uniF,       dis4uniB,       nmat,   n2mat,  &
        nftmx,  nonmx,  nt,     TPV = -1,       output_plastic,     &
        nnx,    nnz,    insertFaultType,    writeCompTime = 0,   outputGroundMotion,&
        surface_nnode = 0,  &
        ! Thickness (counted by nodes) of PML.nPML=6
        nPML = 6, &
        totalNumOfOffSt,  numOfOnFaultStCount,  numOfOffFaultStCount, &
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
     
    real (kind = dp), allocatable, dimension(:) :: nodalForceArr,   nodalMassArr,&
        v1,     stressArr,     miuonf, vponf,  eleporep,       &
        pstrain,eledet, fnms,   fxmin,  fxmax,  fymin,  fymax,  &
        fzmin,  fzmax
    real (kind = dp), allocatable, dimension(:,:) :: meshCoor, &
        dispArr,      velArr,      mat,    localShapeFunc,    fnft,   arn, &
        slp4fri, elemass,ss,     plane1, plane2, &
        Tatnode,patnode,OffFaultStGramSCEC,   material,       rough_geo,      &
        x4nds
    real (kind = dp), allocatable, dimension(:,:,:) :: fric,    &
        un,     us,     ud,     onFaultQuantHistSCECForm, eleshp, phi,    &
        fltxyz, xonfs
    real (kind = dp), allocatable, dimension(:,:,:,:):: onFaultTPHist

    integer (kind = 4), allocatable, dimension(:) :: nftnd,     &
        eqNumIndexArr,    stressCompIndexArr,    elemTypeArr,     eqNumStartIndexLoc,  numOfDofPerNodeArr,   surfaceNodeIdArr,&
        nonfs,  n4yn,   fltl,   fltr,   fltf,   fltb,   fltd,   &
        fltu,   fltgm
    integer (kind = 4), allocatable, dimension(:,:) :: nodeElemIdRelation,     &
        anonfs, idhist, OffFaultStNodeIdIndex
    integer (kind = 4), allocatable, dimension(:,:,:) :: nsmp

end MODULE globalvar           

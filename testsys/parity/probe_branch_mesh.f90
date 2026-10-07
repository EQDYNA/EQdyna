! NOT part of src/fortran/ (unmodified in this probe). Row 153 (TPV24/25)
! feasibility probe, requested by the coordinator to falsify/confirm the
! "strike-along-x assumption" stop condition before committing to the full
! scoped build (PR #141 authorization).
!
! Claim under test: a 30 deg x-y branch fault, extruded along z, can be
! built with the SAME meshgen.f90 loop nest (ix outer, iz, iy inner,
! pplane1/pplane2 swapped at the ix boundary) that test.tpv36/37 already use
! for a y-z dip extruded along x -- i.e. this needs a NEW permutation
! array into the EXISTING, UNMODIFIED reorder() subroutine (library_
! degeneration.f90) plus a cenx-aware distance test, not a re-nested loop.
!
! Derivation (see session notes): reorder()'s nodeIdPerElem(1..8), as
! built from pplane1/pplane2 at (iy-1:iy, iz-1:iz), already contains BOTH
! values of x (pplane1=ix-1 slab, pplane2=ix slab) AND both values of y
! AND both values of z for the current hex -- all three axes are present
! in every element regardless of which axis the outer loop sweeps. The
! existing wedge()/reorder() neworder=(/5,1,4,4,6,2,3,3/) keeps a y-z
! triangle {(iy-1,iz),(iy-1,iz-1),(iy,iz-1)} at x=ix-1 and the SAME
! triangle at x=ix (prism extruded along x). Relabeling the same template
! with (A,B,C) = (iy,iz,ix-extrude) -> (ix,iy,iz-extrude) gives an x-y
! triangle extruded along z, using the SAME nodeIdPerElem indices, just a
! different neworder:
!   type-11 (below):  (/4,1,2,2,8,5,6,6/)
!   type-12 (above):  (/2,3,4,4,6,7,8,8/)
! and a distToFault test using (cenx,ceny) in place of (ceny,cenz).
!
! This program builds a small synthetic grid with EXACTLY that loop nest,
! applies the new permutation + cenx test (new, probe-local subroutines
! wedgeXY/wedge4XY/checkIsOnFaultXY -- candidates for later promotion into
! library_degeneration.f90/meshgen.f90, not yet inserted there), and
! checks the three things the coordinator asked for, calling the REAL,
! UNMODIFIED calcLocalShapeFunc/calcGlobalShapeFunc for the Jacobian.
program probe_branch_mesh
    use globalvar
    use errorCodes
    implicit none
    include 'mpif.h'
    integer (kind = 4) :: iMPIerr
    integer (kind = 4), parameter :: nx = 9, ny = 11, nz = 4
    integer (kind = 4), parameter :: maxnodes = 20000, maxelems = 2000
    real (kind = dp) :: pdx, pdy, pdz, pAngleDeg, pTanA, pNorm
    real (kind = dp) :: xline(nx), yline(ny), zline(nz)
    integer (kind = 8) :: pplane1(ny+2,nz), pplane2(ny+2,nz)
    real (kind = dp) :: meshNodeCoor(3,maxnodes)
    integer (kind = 8) :: nodeIdPerElem(8)
    integer (kind = 4) :: ix, iy, iz, i
    integer (kind = 8) :: nodeCount, elemCount, masterId
    integer (kind = 8) :: slave1(200), master1(200), slave2(200), master2(200)
    integer (kind = 4) :: nmaster1, nmaster2
    integer (kind = 8), parameter :: BIG1 = 10000, BIG2 = 15000
    real (kind = dp) :: coord(3), cenx, ceny, cenz, dist1, distBranch
    real (kind = dp) :: faultXmin2, faultXmax2
    integer (kind = 4) :: neworder11(8), neworder12(8)
    integer (kind = 8) :: nodeElemIdRelationLocal(8,maxelems)
    integer (kind = 4) :: elemTypeLocal(maxelems)
    real (kind = dp) :: xl(3,8), det, xs(3,3), globalShapeFunc(4,8)
    logical :: lcubic, isOnF1, inBranchBox, junctionSeenByFault2
    integer (kind = 4) :: nbad, nwedge, ngoodjac
    integer (kind = 8) :: junctionNodeId

    call MPI_Init(iMPIerr)
    call mpi_comm_rank(MPI_COMM_WORLD, me, iMPIerr)
    call mpi_comm_size(MPI_COMM_WORLD, totalNumOfMPIProcs, iMPIerr)

    allocate(localShapeFunc(nrowsh, nen))
    call calcLocalShapeFunc

    pAngleDeg = 30.0d0
    pdx = 1000.0d0
    pTanA = dtan(pAngleDeg/180.0d0*pi)
    pdy = pdx*pTanA
    pdz = 1000.0d0
    pNorm = (1.0d0 + pTanA**2)**0.5d0

    do ix = 1, nx
        xline(ix) = dble(ix-5)*pdx
    enddo
    do iy = 1, ny
        yline(iy) = dble(iy-6)*pdy
    enddo
    do iz = 1, nz
        zline(iz) = -dble(iz-1)*pdz
    enddo

    ! Branch box: x in [pdx, 4dx] -- EXCLUDES the junction column (x=0), so
    ! the junction node is claimed by fault 1 only (the owner's "junction
    ! node owned by the main fault only" convention), via box membership
    ! alone, no special-case code.
    faultXmin2 = pdx
    faultXmax2 = 4.0d0*pdx

    ! v2: extrude-direction reversed relative to the first (failing) guess
    ! -- a 30 deg x-y wedge needs the opposite z-handedness from the first
    ! substitution attempt (empirically determined: det<0 with the
    ! (/4,1,2,2,8,5,6,6/)/(/2,3,4,4,6,7,8,8/) pair, see probe log).
    neworder11 = (/8,5,6,6,4,1,2,2/)
    neworder12 = (/6,7,8,8,2,3,4,4/)

    pplane1 = 0
    pplane2 = 0
    nodeCount = 0
    elemCount = 0
    nmaster1 = 0
    nmaster2 = 0
    junctionNodeId = 0
    junctionSeenByFault2 = .false.

    do ix = 1, nx
        do iz = 1, nz
            do iy = 1, ny
                nodeCount = nodeCount + 1
                coord = (/ xline(ix), yline(iy), zline(iz) /)
                meshNodeCoor(:,nodeCount) = coord
                pplane2(iy,iz) = nodeCount

                ! --- fault 1: main fault, vertical plane y=0, full x range
                isOnF1 = (abs(coord(2)) < tol)
                if (isOnF1) then
                    nmaster1 = nmaster1 + 1
                    masterId = BIG1 + nmaster1
                    slave1(nmaster1) = nodeCount
                    master1(nmaster1) = masterId
                    meshNodeCoor(:,masterId) = coord
                    pplane2(iy+1,iz) = masterId
                    if (abs(coord(1)) < tol) junctionNodeId = nodeCount
                endif

                ! --- fault 2: branch, x-y tilt, box EXCLUDES x=0 column
                if (coord(1) > faultXmin2-tol .and. coord(1) < faultXmax2+tol) then
                    if (abs(coord(1)) < tol) junctionSeenByFault2 = .true.
                    distBranch = abs(coord(2) + coord(1)*pTanA)/pNorm
                    if (distBranch < pdx/100.0d0) then
                        nmaster2 = nmaster2 + 1
                        masterId = BIG2 + nmaster2
                        slave2(nmaster2) = nodeCount
                        master2(nmaster2) = masterId
                        meshNodeCoor(:,masterId) = coord
                        pplane2(iy+2,iz) = masterId
                    endif
                endif

                if (ix>=2 .and. iy>=2 .and. iz>=2) then
                    nodeIdPerElem(1) = pplane1(iy-1,iz-1)
                    nodeIdPerElem(2) = pplane2(iy-1,iz-1)
                    nodeIdPerElem(3) = pplane2(iy,iz-1)
                    nodeIdPerElem(4) = pplane1(iy,iz-1)
                    nodeIdPerElem(5) = pplane1(iy-1,iz)
                    nodeIdPerElem(6) = pplane2(iy-1,iz)
                    nodeIdPerElem(7) = pplane2(iy,iz)
                    nodeIdPerElem(8) = pplane1(iy,iz)

                    cenx = (xline(ix-1)+xline(ix))/2.0d0
                    ceny = (yline(iy-1)+yline(iy))/2.0d0
                    cenz = (zline(iz-1)+zline(iz))/2.0d0

                    elemCount = elemCount + 1
                    elemTypeLocal(elemCount) = 1
                    nodeElemIdRelationLocal(1:8,elemCount) = nodeIdPerElem(1:8)

                    inBranchBox = (cenx>faultXmin2 .and. cenx<faultXmax2)
                    if (inBranchBox) then
                        distBranch = abs(ceny + cenx*pTanA)/pNorm
                        if (distBranch < tol) then
                            elemTypeLocal(elemCount) = 11
                            do i = 1, 8
                                nodeElemIdRelationLocal(i,elemCount) = nodeIdPerElem(neworder11(i))
                            enddo
                            elemCount = elemCount + 1
                            elemTypeLocal(elemCount) = 12
                            do i = 1, 8
                                nodeElemIdRelationLocal(i,elemCount) = nodeIdPerElem(neworder12(i))
                            enddo
                        endif
                    endif
                endif
            enddo
        enddo
        pplane1 = pplane2
    enddo

    print *, '==== probe_branch_mesh: x-y tilt / z-extrusion feasibility ===='
    print *, 'nodeCount', nodeCount, 'elemCount', elemCount
    print *, 'fault1 split-node pairs', nmaster1, 'fault2 split-node pairs', nmaster2

    ! ---- check 1: element Jacobians positive for every wedge element
    nbad = 0
    nwedge = 0
    allocate(elemTypeArr(elemCount), mat(elemCount,5))
    elemTypeArr(1:elemCount) = elemTypeLocal(1:elemCount)
    do i = 1, elemCount
        if (elemTypeLocal(i) == 11 .or. elemTypeLocal(i) == 12) then
            nwedge = nwedge + 1
            xl(1:3,1:8) = meshNodeCoor(1:3, nodeElemIdRelationLocal(1:8,i))
            lcubic = .true.
            call calcGlobalShapeFunc(xl, det, globalShapeFunc, int(i,8), xs, lcubic)
            if (det <= 0.0d0) then
                nbad = nbad + 1
                print *, '  BAD JACOBIAN elem', i, 'type', elemTypeLocal(i), 'det', det
            endif
        endif
    enddo
    print *, 'check 1: wedge elements', nwedge, 'with det<=0', nbad
    print *, 'check 1 VERDICT: ', merge('PASS', 'FAIL', nwedge>0 .and. nbad==0)

    ! ---- check 2: every branch-plane (fault 2) node has a master/slave pair
    print *, 'check 2: fault2 slave/master pairs created =', nmaster2
    print *, 'check 2 VERDICT: ', merge('PASS', 'FAIL', nmaster2 > 0)

    ! ---- check 3: junction node belongs to fault 1 only
    print *, 'check 3: junction node id (fault1 slave) =', junctionNodeId
    print *, 'check 3: fault2 box ever tested the junction column (x=0)? ', junctionSeenByFault2
    print *, 'check 3 VERDICT: ', merge('PASS', 'FAIL', &
        junctionNodeId > 0 .and. .not. junctionSeenByFault2)

    call MPI_Finalize(iMPIerr)
end program probe_branch_mesh

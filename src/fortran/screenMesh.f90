! Copyright (C) 2006 Benchun Duan, Dunyu Liu
! MIT
subroutine screenMesh(nnode, nel, nx, ny, nz, xline, yline, zline, nftnd0)
    ! Post-meshing screen (owner decision 2026-11-23, board row "wedge-mesh
    ! mis-wiring", deliverable 1). Counterpart: src/python/eqdyna/meshscreen.py
    ! (screen_mesh) -- same five checks, same limits, same message text up to
    ! float formatting. Called from meshgen after the fault-area accumulation
    ! and BEFORE MPI4arn, so arn is this rank's partial sum. A failing mesh
    ! STOPS the run (abortRun, codes 81-85); no switch, no fallback.
    !   sides        a (cell,node) attachment to a pair's MASTER node has the cell
    !                centre on the +un side, one to its SLAVE node on the -un side
    !                (signed cosine strictly > SCREEN_SIDE_MIN); cells whose centre
    !                is beyond the fault patch in x or z are not tested.
    !   bridging     no cell holds both nodes of one pair; pair nodes coincide to tol.
    !   fault_area   arn equals a quarter of the vector areas of the meshed facets
    !                (all-master faces, slave faces where no master side exists)
    !                to SCREEN_AREA_REL; every pair has 1, 2 or 4 facets.
    !   volume_mass  Gauss-point Jacobians > SCREEN_DET_MIN of the cell mean, cell
    !                volumes sum to the rank box volume, interior fault-node volume
    !                shares positive with master/slave ratio <= SCREEN_VOL_RATIO_MAX.
    !   normals      un/us/ud orthonormal; un parallel to the meshed facet normal.
    use globalvar
    use errorCodes
    implicit none
    integer (kind = 8), intent(in) :: nnode, nel
    integer (kind = 4), intent(in) :: nx, ny, nz
    real (kind = dp), intent(in) :: xline(*), yline(*), zline(*)
    integer (kind = 4), intent(in) :: nftnd0(ntotft)

    real (kind = dp), parameter :: SCREEN_SIDE_MIN = 0.0d0
    real (kind = dp), parameter :: SCREEN_AREA_REL = 1.0d-12
    real (kind = dp), parameter :: SCREEN_VOL_REL = 1.0d-10
    real (kind = dp), parameter :: SCREEN_DET_MIN = 0.1d0
    real (kind = dp), parameter :: SCREEN_VOL_RATIO_MAX = 2.0d0
    real (kind = dp), parameter :: SCREEN_NORMAL_TOL = 1.0d-9
    real (kind = dp), parameter :: SCREEN_NORMAL_TOL_ROUGH = 0.2d0
    integer, parameter :: HF(4,6) = reshape((/1,2,3,4, 5,6,7,8, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8/), (/4,6/))
    integer, parameter :: HA(3,8) = reshape((/-1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1, &
                                              -1,-1,1, 1,-1,1, 1,1,1, -1,1,1/), (/3,8/))
    character (len = 12), parameter :: CNAME(5) = (/'sides       ', 'bridging    ', 'fault_area  ', 'volume_mass ', 'normals     '/)
    integer, parameter :: CCODE(5) = (/ERR_MESH_SCREEN_SIDES, ERR_MESH_SCREEN_BRIDGING, ERR_MESH_SCREEN_FAULT_AREA, &
                                       ERR_MESH_SCREEN_VOLUME_MASS, ERR_MESH_SCREEN_NORMALS/)

    integer (kind = 1), allocatable :: role(:)
    integer (kind = 8), allocatable :: pid(:), fcnt(:), c1(:), c2(:)
    integer (kind = 4), allocatable :: pflt(:), pidx(:), poff(:)
    logical, allocatable :: onEdge(:)
    real (kind = dp), allocatable :: share(:), e1(:), e2(:), n1v(:,:), n2v(:,:), fexp(:), relA(:), dev(:), ang(:), &
                                     nacc(:,:), lrat(:), ratio(:), vs(:), vm(:)
    real (kind = dp) :: xlo(ntotft), xhi(ntotft), zlo(ntotft), zhi(ntotft)
    character (len = 4000) :: msg(5), meas
    logical :: failed(5)
    integer (kind = 8) :: e, nft, p, nd, c(8), fn(4), nbadj, nwrong, natt, nmw, nsw, nbr, firstBr, firstBadJ
    integer (kind = 8) :: wE, wN, wP, nOut, nEmpty, firstEmpty, nBadCnt, nNoFacet, firstNoFacet, firstBadCnt
    integer (kind = 4) :: ift, i, j, k, f, rl, a, b, g1, g2, g3, ix, iz, wI, wF
    real (kind = dp) :: cen(3), dv(3), un3(3), va(3), ar, sg, margin, bestM, minM, dn, nv, cosd, vol, dmin, dj(8), jac(3,3), dnn(8,3), &
                        ff(8,3), s3(3), g, relDet, minRelDet, sep, sepMax, relMax, devMax, angMax, vsum, vexp, hh(nx,nz), hc, &
                        lim, lrMax, volRel, nn, cr(3), nf(3), t, maxRatio
    integer :: nEmptyI
    logical :: hasTouch, inside

    g = 1.0d0/sqrt(3.0d0)
    nft = sum(int(nftnd0, 8))
    failed = .false.
    msg = ' '
    meas = ' '
    allocate(role(nnode), pid(nnode), share(nnode))
    role = 0; pid = 0; share = 0.0d0
    allocate(pflt(max(nft,1_8)), pidx(max(nft,1_8)), poff(ntotft), onEdge(max(nft,1_8)))
    poff = 0
    do ift = 2, ntotft
        poff(ift) = poff(ift-1) + nftnd0(ift-1)
    enddo
    xlo = huge(1.0d0); xhi = -huge(1.0d0); zlo = huge(1.0d0); zhi = -huge(1.0d0)
    do ift = 1, ntotft
        do i = 1, nftnd0(ift)
            p = poff(ift) + i
            role(nsmp(1,i,ift)) = 1
            role(nsmp(2,i,ift)) = 2
            pid(nsmp(1,i,ift)) = p
            pid(nsmp(2,i,ift)) = p
            pflt(p) = ift
            pidx(p) = i
            xlo(ift) = min(xlo(ift), meshCoor(1,nsmp(2,i,ift)))
            xhi(ift) = max(xhi(ift), meshCoor(1,nsmp(2,i,ift)))
            zlo(ift) = min(zlo(ift), meshCoor(3,nsmp(2,i,ift)))
            zhi(ift) = max(zhi(ift), meshCoor(3,nsmp(2,i,ift)))
        enddo
        do i = 1, nftnd0(ift)
            p = poff(ift) + i
            nd = nsmp(2,i,ift)
            onEdge(p) = abs(meshCoor(1,nd)-xlo(ift)) <= tol .or. abs(meshCoor(1,nd)-xhi(ift)) <= tol .or. &
                        abs(meshCoor(3,nd)-zlo(ift)) <= tol .or. abs(meshCoor(3,nd)-zhi(ift)) <= tol
        enddo
    enddo
    allocate(fcnt(max(nft,1_8)), c1(max(nft,1_8)), c2(max(nft,1_8)), e1(max(nft,1_8)), e2(max(nft,1_8)), &
             n1v(3,max(nft,1_8)), n2v(3,max(nft,1_8)), fexp(max(nft,1_8)), relA(max(nft,1_8)), dev(max(nft,1_8)), &
             ang(max(nft,1_8)), nacc(3,max(nft,1_8)), lrat(max(nft,1_8)), ratio(max(nft,1_8)), vs(max(nft,1_8)), vm(max(nft,1_8)))
    c1 = 0; c2 = 0; e1 = 0.0d0; e2 = 0.0d0; n1v = 0.0d0; n2v = 0.0d0

    natt = 0; nwrong = 0; nmw = 0; nsw = 0; minM = huge(1.0d0); bestM = huge(1.0d0); wE = 0; wN = 0; wP = 0
    nbr = 0; firstBr = 0; nbadj = 0; firstBadJ = 0; minRelDet = huge(1.0d0); vsum = 0.0d0

    do e = 1, nel
        c = nodeElemIdRelation(1:8,e)
        ! ---- volume_mass: Gauss-point Jacobians, volume, node shares (all cells)
        k = 0
        do g3 = -1, 1, 2
            do g2 = -1, 1, 2
                do g1 = -1, 1, 2
                    k = k + 1
                    s3 = (/g1*g, g2*g, g3*g/)
                    do i = 1, 8
                        do j = 1, 3
                            ff(i,j) = 1.0d0 + HA(j,i)*s3(j)
                        enddo
                        dnn(i,1) = HA(1,i)*ff(i,2)*ff(i,3)/8.0d0
                        dnn(i,2) = HA(2,i)*ff(i,1)*ff(i,3)/8.0d0
                        dnn(i,3) = HA(3,i)*ff(i,1)*ff(i,2)/8.0d0
                    enddo
                    do a = 1, 3
                        do b = 1, 3
                            jac(a,b) = 0.0d0
                            do i = 1, 8
                                jac(a,b) = jac(a,b) + dnn(i,a)*meshCoor(b,c(i))
                            enddo
                        enddo
                    enddo
                    dj(k) = jac(1,1)*(jac(2,2)*jac(3,3)-jac(2,3)*jac(3,2)) - jac(1,2)*(jac(2,1)*jac(3,3)-jac(2,3)*jac(3,1)) &
                          + jac(1,3)*(jac(2,1)*jac(3,2)-jac(2,2)*jac(3,1))
                enddo
            enddo
        enddo
        vol = sum(dj)
        dmin = minval(dj)
        if (vol > 0.0d0) then
            relDet = dmin/(vol/8.0d0)
        else
            relDet = dmin
        endif
        minRelDet = min(minRelDet, relDet)
        if (.not. (relDet > SCREEN_DET_MIN)) then
            nbadj = nbadj + 1
            if (nbadj == 1) firstBadJ = e
        endif
        vsum = vsum + vol
        do i = 1, 8
            share(c(i)) = share(c(i)) + vol/8.0d0
        enddo

        hasTouch = .false.
        do i = 1, 8
            if (role(c(i)) /= 0) hasTouch = .true.
        enddo
        if (.not. hasTouch) cycle
        do j = 1, 3
            cen(j) = 0.0d0
            do i = 1, 8
                cen(j) = cen(j) + meshCoor(j,c(i))
            enddo
            cen(j) = cen(j)/8.0d0
        enddo
        ! ---- sides
        do i = 1, 8
            if (role(c(i)) == 0) cycle
            if (any(c(1:i-1) == c(i))) cycle
            p = pid(c(i)); ift = pflt(p)
            inside = cen(1) >= xlo(ift) .and. cen(1) <= xhi(ift) .and. cen(3) >= zlo(ift) .and. cen(3) <= zhi(ift)
            if (.not. inside) cycle
            dv = cen - meshCoor(1:3,c(i))
            dn = dot_product(dv, un(1:3,pidx(p),ift))
            nv = sqrt(dot_product(dv, dv))
            if (nv > 0.0d0) then
                cosd = dn/nv
            else
                cosd = dn
            endif
            if (role(c(i)) == 2) then
                margin = cosd
            else
                margin = -cosd
            endif
            natt = natt + 1
            minM = min(minM, margin)
            if (.not. (margin > SCREEN_SIDE_MIN)) then
                nwrong = nwrong + 1
                if (role(c(i)) == 2) then
                    nmw = nmw + 1
                else
                    nsw = nsw + 1
                endif
                if (margin < bestM .or. (margin == bestM .and. e == wE .and. c(i) < wN)) then
                    bestM = margin; wE = e; wN = c(i); wP = p
                endif
            endif
        enddo
        ! ---- bridging
        do i = 1, 8
            if (role(c(i)) /= 1) cycle
            do j = 1, 8
                if (role(c(j)) == 2 .and. pid(c(j)) == pid(c(i))) then
                    nbr = nbr + 1
                    if (nbr == 1) firstBr = e
                    goto 100
                endif
            enddo
        enddo
100     continue
        ! ---- facets
        do f = 1, 6
            fn = c(HF(1:4,f))
            if (fn(1)==fn(2) .or. fn(1)==fn(3) .or. fn(1)==fn(4) .or. fn(2)==fn(3) .or. fn(2)==fn(4) .or. fn(3)==fn(4)) cycle
            do rl = 2, 1, -1
                if (.not. (role(fn(1))==rl .and. role(fn(2))==rl .and. role(fn(3))==rl .and. role(fn(4))==rl)) cycle
                va = 0.5d0*cross3(meshCoor(1:3,fn(3))-meshCoor(1:3,fn(1)), meshCoor(1:3,fn(4))-meshCoor(1:3,fn(2)))
                ar = sqrt(dot_product(va, va))
                do k = 1, 4
                    p = pid(fn(k))
                    sg = dot_product(va, un(1:3,pidx(p),pflt(p)))
                    if (sg > 0.0d0) then
                        sg = 1.0d0
                    elseif (sg < 0.0d0) then
                        sg = -1.0d0
                    else
                        sg = 1.0d0
                    endif
                    if (rl == 2) then
                        c2(p) = c2(p) + 1; e2(p) = e2(p) + 0.25d0*ar; n2v(1:3,p) = n2v(1:3,p) + va*sg
                    else
                        c1(p) = c1(p) + 1; e1(p) = e1(p) + 0.25d0*ar; n1v(1:3,p) = n1v(1:3,p) + va*sg
                    endif
                enddo
            enddo
        enddo
    enddo

    ! ---- sides message
    if (nwrong > 0) then
        failed(1) = .true.
        msg(1) = trim(ii(nwrong))//' of '//trim(ii(natt))//' (cell,node) attachments on the wrong side of the fault plane ('// &
            trim(ii(nmw))//' master-attached on the footwall side, '//trim(ii(nsw))// &
            ' slave-attached on the hanging-wall side); worst: element '//trim(ii(wE))//' (type '// &
            trim(ii(int(elemTypeArr(wE),8)))//') holds '//trim(merge('master', 'slave ', role(wN)==2))//' node '//trim(ii(wN))// &
            ' of pair '//trim(ii(int(pidx(wP),8)))//' (fault '//trim(ii(int(pflt(wP),8)))//'), margin '//trim(fe(bestM))
    endif
    ! ---- bridging
    sepMax = -1.0d0; wI = 0; wF = 1
    do ift = 1, ntotft
        do i = 1, nftnd0(ift)
            sep = maxval(abs(meshCoor(1:3,nsmp(1,i,ift)) - meshCoor(1:3,nsmp(2,i,ift))))
            if (sep > sepMax) then
                sepMax = sep; wI = i; wF = ift
            endif
        enddo
    enddo
    if (nft == 0) sepMax = 0.0d0
    if (nbr > 0) then
        failed(2) = .true.
        msg(2) = trim(ii(nbr))//' cell(s) hold both nodes of a split pair, first element '//trim(ii(firstBr))//' (type '// &
                 trim(ii(int(elemTypeArr(firstBr),8)))//')'
    endif
    if (nft > 0) then
        if (.not. (sepMax <= tol)) then
            failed(2) = .true.
            if (len_trim(msg(2)) > 0) msg(2) = trim(msg(2))//' ; '
            msg(2) = trim(msg(2))//'pair '//trim(ii(int(poff(wF)+wI,8)))//' (fault '//trim(ii(int(wF,8)))//') nodes '// &
                trim(ii(nsmp(1,wI,wF)))//'/'//trim(ii(nsmp(2,wI,wF)))//' are '//trim(fe(sepMax))//' m apart (tol '//trim(fe(tol))//')'
        endif
    endif
    ! ---- facet choice, fault_area
    do p = 1, nft
        if (c2(p) > 0) then
            fcnt(p) = c2(p); fexp(p) = e2(p); nacc(1:3,p) = n2v(1:3,p)
        else
            fcnt(p) = c1(p); fexp(p) = e1(p); nacc(1:3,p) = n1v(1:3,p)
        endif
    enddo
    nNoFacet = 0; firstNoFacet = 0; nBadCnt = 0; firstBadCnt = 0
    relMax = 0.0d0
    do p = 1, nft
        if (fcnt(p) == 0) then
            nNoFacet = nNoFacet + 1
            if (nNoFacet == 1) firstNoFacet = p
        elseif (fcnt(p) /= 1 .and. fcnt(p) /= 2 .and. fcnt(p) /= 4) then
            nBadCnt = nBadCnt + 1
            if (nBadCnt == 1) firstBadCnt = p
        endif
        if (fcnt(p) > 0) then
            if (fexp(p) > 0.0d0) then
                relA(p) = abs(arn(pidx(p),pflt(p)) - fexp(p))/fexp(p)
            else
                relA(p) = abs(arn(pidx(p),pflt(p)) - fexp(p))
            endif
        else
            relA(p) = 0.0d0
        endif
        relMax = max(relMax, relA(p))
    enddo
    if (nNoFacet > 0) then
        failed(3) = .true.
        msg(3) = trim(ii(nNoFacet))//' pair(s) have no meshed facet on either side, first pair '// &
                 trim(ii(int(pidx(firstNoFacet),8)))//' (fault '//trim(ii(int(pflt(firstNoFacet),8)))//')'
    endif
    if (nBadCnt > 0) then
        failed(3) = .true.
        if (len_trim(msg(3)) > 0) msg(3) = trim(msg(3))//' ; '
        msg(3) = trim(msg(3))//trim(ii(nBadCnt))//' pair(s) carry a facet count other than 1, 2 or 4, first pair '// &
                 trim(ii(int(pidx(firstBadCnt),8)))//' has '//trim(ii(fcnt(firstBadCnt)))
    endif
    do ift = 1, ntotft
        if (nftnd0(ift) == 0) cycle
        wP = poff(ift) + 1; nOut = 0
        do p = poff(ift)+1, poff(ift)+nftnd0(ift)
            if (relA(p) > relA(wP)) wP = p
            if (relA(p) > SCREEN_AREA_REL) nOut = nOut + 1
        enddo
        if (.not. (relA(wP) <= SCREEN_AREA_REL)) then
            failed(3) = .true.
            if (len_trim(msg(3)) > 0) msg(3) = trim(msg(3))//' ; '
            msg(3) = trim(msg(3))//'fault '//trim(ii(int(ift,8)))//': node area of pair '//trim(ii(int(pidx(wP),8)))//' is '// &
                trim(fe(arn(pidx(wP),ift)))//' but the incident facets give '//trim(fe(fexp(wP)))//' (rel '//trim(fe(relA(wP)))// &
                ' > '//trim(fe(SCREEN_AREA_REL))//'); '//trim(ii(nOut))//' of '//trim(ii(int(nftnd0(ift),8)))//' pairs off'
        endif
    enddo
    ! ---- normals
    devMax = 0.0d0; angMax = 0.0d0
    do p = 1, nft
        ift = pflt(p); i = pidx(p)
        dev(p) = max(abs(sqrt(dot_product(un(1:3,i,ift),un(1:3,i,ift)))-1.0d0), abs(sqrt(dot_product(us(1:3,i,ift),us(1:3,i,ift)))-1.0d0), &
                     abs(sqrt(dot_product(ud(1:3,i,ift),ud(1:3,i,ift)))-1.0d0), abs(dot_product(un(1:3,i,ift),us(1:3,i,ift))), &
                     abs(dot_product(un(1:3,i,ift),ud(1:3,i,ift))), abs(dot_product(us(1:3,i,ift),ud(1:3,i,ift))))
        nn = sqrt(dot_product(nacc(1:3,p), nacc(1:3,p)))
        if (nn > 0.0d0) then
            nf = nacc(1:3,p)/nn
        else
            nf = nacc(1:3,p)
        endif
        cr = cross3(un(1:3,i,ift), nf)
        if (fcnt(p) > 0) then
            ang(p) = sqrt(dot_product(cr, cr))
        else
            ang(p) = 0.0d0
        endif
        devMax = max(devMax, dev(p)); angMax = max(angMax, ang(p))
    enddo
    do ift = 1, ntotft
        if (nftnd0(ift) == 0) cycle
        if (insertFaultType == 0) then
            lim = SCREEN_NORMAL_TOL
        else
            lim = SCREEN_NORMAL_TOL_ROUGH
        endif
        wP = poff(ift) + 1
        do p = poff(ift)+1, poff(ift)+nftnd0(ift)
            if (dev(p) > dev(wP)) wP = p
        enddo
        if (.not. (dev(wP) <= SCREEN_NORMAL_TOL)) then
            failed(5) = .true.
            if (len_trim(msg(5)) > 0) msg(5) = trim(msg(5))//' ; '
            msg(5) = trim(msg(5))//'fault '//trim(ii(int(ift,8)))//': (un,us,ud) of pair '//trim(ii(int(pidx(wP),8)))// &
                ' is not orthonormal (worst deviation '//trim(fe(dev(wP)))//' > '//trim(fe(SCREEN_NORMAL_TOL))//')'
        endif
        wP = poff(ift) + 1; nOut = 0
        do p = poff(ift)+1, poff(ift)+nftnd0(ift)
            if (ang(p) > ang(wP)) wP = p
            if (ang(p) > lim) nOut = nOut + 1
        enddo
        if (.not. (ang(wP) <= lim)) then
            failed(5) = .true.
            if (len_trim(msg(5)) > 0) msg(5) = trim(msg(5))//' ; '
            msg(5) = trim(msg(5))//'fault '//trim(ii(int(ift,8)))//': un of pair '//trim(ii(int(pidx(wP),8)))// &
                ' is not perpendicular to its meshed facets (sin angle '//trim(fe(ang(wP)))//' > '//trim(fe(lim))//'); '// &
                trim(ii(nOut))//' of '//trim(ii(int(nftnd0(ift),8)))//' pairs off'
        endif
    enddo
    ! ---- volume_mass
    if (nbadj > 0) then
        failed(4) = .true.
        msg(4) = trim(ii(nbadj))//' of '//trim(ii(nel))//' cells have a Gauss-point Jacobian below '//trim(fe(SCREEN_DET_MIN))// &
                 ' of their mean, first element '//trim(ii(firstBadJ))//' (type '//trim(ii(int(elemTypeArr(firstBadJ),8)))// &
                 ', volume '//trim(fe(vol_of(firstBadJ)))//')'
    endif
    do ix = 1, nx
        do iz = 1, nz
            hh(ix,iz) = meshCoor(2,((ix-1)*nz+(iz-1))*ny+ny) - meshCoor(2,((ix-1)*nz+(iz-1))*ny+1)
        enddo
    enddo
    vexp = 0.0d0
    do ix = 1, nx-1
        do iz = 1, nz-1
            hc = 0.25d0*(hh(ix,iz) + hh(ix+1,iz) + hh(ix,iz+1) + hh(ix+1,iz+1))
            vexp = vexp + (xline(ix+1)-xline(ix))*(zline(iz+1)-zline(iz))*hc
        enddo
    enddo
    volRel = abs(vsum - vexp)/vexp
    if (.not. (volRel <= SCREEN_VOL_REL)) then
        failed(4) = .true.
        if (len_trim(msg(4)) > 0) msg(4) = trim(msg(4))//' ; '
        msg(4) = trim(msg(4))//'cell volumes sum to '//trim(fe(vsum))//' m^3 but the rank box holds '//trim(fe(vexp))// &
                 ' m^3 (rel '//trim(fe(volRel))//' > '//trim(fe(SCREEN_VOL_REL))//')'
    endif
    nEmpty = 0; firstEmpty = 0; nOut = 0; lrMax = 0.0d0; wP = 1
    do ift = 1, ntotft
        do i = 1, nftnd0(ift)
            p = poff(ift) + i
            vs(p) = share(nsmp(1,i,ift)); vm(p) = share(nsmp(2,i,ift))
            j = int(mod(nsmp(1,i,ift)-1, int(ny,8)))
            if ((vs(p) <= 0.0d0 .or. vm(p) <= 0.0d0) .and. .not. (j == 0 .or. j == ny-1)) then
                nEmpty = nEmpty + 1
                if (nEmpty == 1) firstEmpty = p
            endif
            if (vs(p) > 0.0d0 .and. vm(p) > 0.0d0 .and. .not. onEdge(p)) then
                ratio(p) = vm(p)/vs(p)
            else
                ratio(p) = 1.0d0
            endif
            t = abs(log(ratio(p)))
            if (t > lrMax) then
                lrMax = t; wP = p
            endif
            if (t > log(SCREEN_VOL_RATIO_MAX)) nOut = nOut + 1
        enddo
    enddo
    if (nEmpty > 0) then
        failed(4) = .true.
        if (len_trim(msg(4)) > 0) msg(4) = trim(msg(4))//' ; '
        msg(4) = trim(msg(4))//trim(ii(nEmpty))//' fault node(s) have no volume share on a side, first pair '// &
                 trim(ii(int(pidx(firstEmpty),8)))//' (fault '//trim(ii(int(pflt(firstEmpty),8)))//')'
    endif
    if (nOut > 0) then
        failed(4) = .true.
        if (len_trim(msg(4)) > 0) msg(4) = trim(msg(4))//' ; '
        msg(4) = trim(msg(4))//'fault node pair '//trim(ii(int(pidx(wP),8)))//' (fault '//trim(ii(int(pflt(wP),8)))// &
                 ') has master/slave volume share '//trim(fe(ratio(wP)))//' (limit '//trim(fe(SCREEN_VOL_RATIO_MAX))//'); '// &
                 trim(ii(nOut))//' pairs outside'
    endif
    maxRatio = exp(lrMax)

    meas = 'sides_attachments='//trim(ii(natt))//' sides_min_margin='//trim(fe(minM))//' sides_wrong_master='//trim(ii(nmw))// &
           ' sides_wrong_slave='//trim(ii(nsw))//' bridging_cells='//trim(ii(nbr))//' coincidence_max='//trim(fe(sepMax))// &
           ' area_nofacet_pairs='//trim(ii(nNoFacet))//' area_badcount_pairs='//trim(ii(nBadCnt))//' area_facet_rel_max='// &
           trim(fe(relMax))//' normal_orthonormal_dev_max='//trim(fe(devMax))//' normal_facet_sin_max='//trim(fe(angMax))// &
           ' det_min_rel='//trim(fe(minRelDet))//' det_bad_cells='//trim(ii(nbadj))//' volume_rel='//trim(fe(volRel))// &
           ' share_empty_pairs='//trim(ii(nEmpty))//' mass_ratio_max='//trim(fe(maxRatio))//' mass_ratio_out_pairs='//trim(ii(nOut))

    if (any(failed)) then
        block
            character (len = 100) :: names
            character (len = 20000) :: reason
            integer :: code
            names = ' '
            reason = 'mesh screening FAILED: '
            code = 0
            do k = 1, 5
                if (.not. failed(k)) cycle
                if (code == 0) code = CCODE(k)
                if (len_trim(names) > 0) names = trim(names)//','
                names = trim(names)//trim(CNAME(k))
                if (k > 1 .and. len_trim(reason) > 23) reason = trim(reason)//' | '
                reason = trim(reason)//'['//trim(CNAME(k))//'] '//trim(msg(k))
            enddo
            write(*,'(a)') 'mesh-screen: FAIL checks='//trim(names)//' '//trim(meas)
            call abortRun(code, trim(reason))
        end block
    endif
    if (meshOnly) write(*,'(a)') 'mesh-screen: PASS checks=none '//trim(meas)
    deallocate(role, pid, share, pflt, pidx, poff, onEdge, fcnt, c1, c2, e1, e2, n1v, n2v, fexp, relA, dev, ang, nacc, &
               lrat, ratio, vs, vm)

contains
    function cross3(u, v) result(w)
        real (kind = dp), intent(in) :: u(3), v(3)
        real (kind = dp) :: w(3)
        w(1) = u(2)*v(3) - u(3)*v(2)
        w(2) = u(3)*v(1) - u(1)*v(3)
        w(3) = u(1)*v(2) - u(2)*v(1)
    end function cross3
    function ii(n) result(s)
        integer (kind = 8), intent(in) :: n
        character (len = 24) :: s
        write(s,'(I0)') n
    end function ii
    function fe(v) result(s)
        real (kind = dp), intent(in) :: v
        character (len = 14) :: s
        write(s,'(ES11.3)') v
        s = adjustl(s)
    end function fe
    function vol_of(el) result(v)
        integer (kind = 8), intent(in) :: el
        real (kind = dp) :: v
        integer :: q1, q2, q3, kk, i1, j1, a1, b1
        real (kind = dp) :: s(3), f2(8,3), d2(8,3), jc(3,3)
        integer (kind = 8) :: cc(8)
        cc = nodeElemIdRelation(1:8,el)
        v = 0.0d0
        do q3 = -1, 1, 2
            do q2 = -1, 1, 2
                do q1 = -1, 1, 2
                    s = (/q1*g, q2*g, q3*g/)
                    do i1 = 1, 8
                        do j1 = 1, 3
                            f2(i1,j1) = 1.0d0 + HA(j1,i1)*s(j1)
                        enddo
                        d2(i1,1) = HA(1,i1)*f2(i1,2)*f2(i1,3)/8.0d0
                        d2(i1,2) = HA(2,i1)*f2(i1,1)*f2(i1,3)/8.0d0
                        d2(i1,3) = HA(3,i1)*f2(i1,1)*f2(i1,2)/8.0d0
                    enddo
                    do a1 = 1, 3
                        do b1 = 1, 3
                            jc(a1,b1) = 0.0d0
                            do i1 = 1, 8
                                jc(a1,b1) = jc(a1,b1) + d2(i1,a1)*meshCoor(b1,cc(i1))
                            enddo
                        enddo
                    enddo
                    v = v + jc(1,1)*(jc(2,2)*jc(3,3)-jc(2,3)*jc(3,2)) - jc(1,2)*(jc(2,1)*jc(3,3)-jc(2,3)*jc(3,1)) &
                          + jc(1,3)*(jc(2,1)*jc(3,2)-jc(2,2)*jc(3,1))
                enddo
            enddo
        enddo
    end function vol_of
end subroutine screenMesh

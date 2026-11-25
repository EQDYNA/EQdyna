! Copyright (C) 2006 Benchun Duan, Dunyu Liu
! MIT
subroutine screenMesh(nnode, nel, nftnd0)
    ! Post-meshing screen (board row 167). Counterpart: src/python/eqdyna/meshscreen.py
    ! (screen_mesh): same two checks, same counts. Called from meshgen after the
    ! fault-area accumulation and BEFORE MPI4arn, on this rank's box. A failing
    ! mesh STOPS the run (abortRun, codes 81-82); no switch, no fallback.
    !   sides     a (cell,node) attachment to a pair's MASTER node has the cell
    !             centre on the +un side, one to its SLAVE node on the -un side
    !             (signed cosine strictly > 0). Beyond the fault patch (cell centre
    !             outside it in x or z) only MASTER attachments are tested; slave
    !             attachments there are exempt (edge-weld convention).
    !   bridging  no cell holds both nodes of one pair; pair nodes coincide to tol.
    use globalvar
    use errorCodes
    implicit none
    integer (kind = 8), intent(in) :: nnode, nel
    integer (kind = 4), intent(in) :: nftnd0(ntotft)

    integer (kind = 1), allocatable :: role(:)
    integer (kind = 4), allocatable :: pn(:), pf(:)
    real (kind = dp) :: xlo(ntotft), xhi(ntotft), zlo(ntotft), zhi(ntotft)
    character (len = 2000) :: msgSides, msgBridge
    integer (kind = 8) :: e, c(8), s, m, nft, natt, nwrong, nmw, nsw, nbr, firstWrong, firstBr, wS, wM
    integer (kind = 4) :: ift, i, j, k
    real (kind = dp) :: cen(3), dv(3), nv, cosd, margin, minM, sep, sepMax
    logical :: br

    nft = sum(int(nftnd0, 8))
    allocate(role(nnode), pn(nnode), pf(nnode))
    role = 0; pn = 0; pf = 0
    xlo = huge(1.0d0); xhi = -huge(1.0d0); zlo = huge(1.0d0); zhi = -huge(1.0d0)
    sepMax = 0.0d0; wS = 0; wM = 0
    do ift = 1, ntotft
        do i = 1, nftnd0(ift)
            s = nsmp(1,i,ift); m = nsmp(2,i,ift)
            role(s) = 1; role(m) = 2
            pn(s) = i; pn(m) = i; pf(s) = ift; pf(m) = ift
            xlo(ift) = min(xlo(ift), meshCoor(1,m)); xhi(ift) = max(xhi(ift), meshCoor(1,m))
            zlo(ift) = min(zlo(ift), meshCoor(3,m)); zhi(ift) = max(zhi(ift), meshCoor(3,m))
            sep = maxval(abs(meshCoor(1:3,s) - meshCoor(1:3,m)))
            if (sep > sepMax) then
                sepMax = sep; wS = s; wM = m
            endif
        enddo
    enddo

    natt = 0; nwrong = 0; nmw = 0; nsw = 0; nbr = 0; firstWrong = 0; firstBr = 0; minM = huge(1.0d0)
    do e = 1, nel
        c = nodeElemIdRelation(1:8,e)
        if (all(role(c) == 0)) cycle
        do j = 1, 3
            cen(j) = 0.0d0
            do i = 1, 8
                cen(j) = cen(j) + meshCoor(j,c(i))
            enddo
            cen(j) = cen(j)/8.0d0
        enddo
        br = .false.
        do i = 1, 8
            if (role(c(i)) == 0) cycle
            if (any(c(1:i-1) == c(i))) cycle
            ift = pf(c(i))
            if (role(c(i)) == 1) then
                do k = 1, 8
                    if (role(c(k)) == 2 .and. pf(c(k)) == ift .and. pn(c(k)) == pn(c(i))) br = .true.
                enddo
            endif
            ! Slave attachments beyond the patch are exempt (fxmax+dx edge-weld convention; same measured
            ! reason as meshscreen.py); master attachments are tested everywhere (row 169 (b)).
            if (role(c(i)) == 1 .and. .not. (cen(1) >= xlo(ift) .and. cen(1) <= xhi(ift) .and. &
                cen(3) >= zlo(ift) .and. cen(3) <= zhi(ift))) cycle
            dv = cen - meshCoor(1:3,c(i))
            nv = sqrt(dot_product(dv, dv))
            cosd = dot_product(dv, un(1:3,pn(c(i)),ift))
            if (nv > 0.0d0) cosd = cosd/nv
            margin = cosd
            if (role(c(i)) == 1) margin = -cosd
            natt = natt + 1
            minM = min(minM, margin)
            if (.not. (margin > 0.0d0)) then
                nwrong = nwrong + 1
                if (role(c(i)) == 2) nmw = nmw + 1
                if (role(c(i)) == 1) nsw = nsw + 1
                if (firstWrong == 0) firstWrong = e
            endif
        enddo
        if (br) then
            nbr = nbr + 1
            if (nbr == 1) firstBr = e
        endif
    enddo

    msgSides = ' '; msgBridge = ' '
    if (nwrong > 0) then
        msgSides = trim(ii(nwrong))//' of '//trim(ii(natt))//' (cell,node) attachments on the wrong side of the fault plane ('// &
            trim(ii(nmw))//' master-attached on the footwall side, '//trim(ii(nsw))//' slave-attached on the hanging-wall side)'// &
            ', first element '//trim(ii(firstWrong))//' (type '//trim(ii(int(elemTypeArr(firstWrong),8)))//')'
    endif
    if (nbr > 0) then
        msgBridge = trim(ii(nbr))//' cell(s) hold both nodes of a split pair, first element '//trim(ii(firstBr))// &
            ' (type '//trim(ii(int(elemTypeArr(firstBr),8)))//')'
    endif
    if (.not. (sepMax <= tol)) then
        if (nbr > 0) msgBridge = trim(msgBridge)//' ; '
        msgBridge = trim(msgBridge)//'nodes '//trim(ii(wS))//'/'//trim(ii(wM))//' of a pair are '//trim(fe(sepMax))// &
            ' m apart (tol '//trim(fe(tol))//')'
    endif

    block
        character (len = 400) :: meas
        character (len = 6000) :: reason
        character (len = 20) :: names
        integer :: code
        meas = 'elements='//trim(ii(nel))//' nodes='//trim(ii(nnode))//' fault_pairs='//trim(ii(nft))// &
               ' sides_attachments='//trim(ii(natt))//' sides_min_margin='//trim(fe(minM))//' sides_wrong_master='// &
               trim(ii(nmw))//' sides_wrong_slave='//trim(ii(nsw))//' bridging_cells='//trim(ii(nbr))// &
               ' coincidence_max='//trim(fe(sepMax))
        names = ' '; reason = ' '; code = 0
        if (len_trim(msgSides) > 0) then
            code = ERR_MESH_SCREEN_SIDES
            names = 'sides'
            reason = '[sides] '//trim(msgSides)
        endif
        if (len_trim(msgBridge) > 0) then
            if (code == 0) then
                code = ERR_MESH_SCREEN_BRIDGING
                names = 'bridging'
                reason = '[bridging] '//trim(msgBridge)
            else
                names = 'sides,bridging'
                reason = trim(reason)//' | [bridging] '//trim(msgBridge)
            endif
        endif
        if (code /= 0) then
            write(*,'(a)') 'mesh-screen: FAIL checks='//trim(names)//' '//trim(meas)
            call abortRun(code, 'mesh screening FAILED: '//trim(reason))
        endif
        ! --mesh-only: every rank prints; a normal run: rank 0 only (row 169 (e)), so a solve shows the screen ran.
        if (meshOnly .or. me == 0) write(*,'(a)') 'mesh-screen: PASS checks=none '//trim(meas)
    end block
    deallocate(role, pn, pf)

contains
    function ii(n) result(str)
        integer (kind = 8), intent(in) :: n
        character (len = 24) :: str
        write(str,'(I0)') n
    end function ii
    function fe(v) result(str)
        real (kind = dp), intent(in) :: v
        character (len = 14) :: str
        write(str,'(ES11.3)') v
        str = adjustl(str)
    end function fe
end subroutine screenMesh

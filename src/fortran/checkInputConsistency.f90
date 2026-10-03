subroutine checkInputConsistency

    use globalvar
    use errorCodes
    implicit none
    integer (kind = 4) :: i, j
    character (len = 12) :: itoa
    character (len = 24) :: dtoa

    if (C_elastic==0.and.C_Q==1) then
        call abortRun(ERR_CFG_Q_NEEDS_ELASTIC, &
            'Q model (C_Q=1) can only work with the elastic code (C_elastic=1). Set C_Q=0 or C_elastic=1.')
    endif
    if (C_Q==1.and.rat>1) then
        call abortRun(ERR_CFG_Q_NEEDS_UNIFORM, &
            'Q model (C_Q=1) can only work with uniform element size; rat must be 1.0.')
    endif
    if (output_plastic == 1 .and. C_elastic/=0) then
        write(*,*) 'Now, C_elastic = ', C_elastic
        call abortRun(ERR_CFG_PLASTIC_OUTPUT, &
            'Plastic strains are only output for C_elastic=0. Set output_plastic=0 or C_elastic=0.')
    endif

    ! Row 17 (multi-fault) guards, modelled on eqquasi's eeac6f9 (fault plane
    ! must land on a mesh node line) and a761f33 (x/z bound commensurability).
    ! No-op at ntotft==1 (every loop below is over a single fault, and the
    ! i<j distinctness loop does not execute for ntotft<2).
    ! These multi-fault guards are about checkIsOnFault's PLANAR,
    ! C_degen==0 branch (one y-plane per fault, tested by exact-ish equality
    ! against fltxyz(1,2,iFault)) -- the narrow two-vertical-parallel-faults
    ! scope this release targets. C_degen>3 (wedge degeneration, e.g.
    ! test.tpv36/test.tpv37) is a DIFFERENT, pre-existing, already
    ! single-fault-only dipping-fault mechanism with a legitimate
    ! fymin /= fymax (the fault spans a range of grid y by construction) --
    ! orthogonal to this work and explicitly out of scope, so these guards
    ! do not apply to it.
    !
    ! Row 17 rebased (restore per-fault mesh extent, invent nothing): the
    ! "y outside the fixed belt" and "x/z must equal fault 1's" checks that
    ! used to live here are GONE, not relaxed -- they were refusing exactly
    ! the two limitations getLocalOneDimCoorArrAndSize (meshgen.f90) no
    ! longer has. The belt is now derived FROM every fault's own box (union
    ! min/max), so a fault's y can never fall outside it, and independent
    ! per-fault x/z extents are now meshed, not refused. What replaces both
    ! checks is one generic commensurability guard inside
    ! getLocalOneDimCoorArrAndSize itself (it alone knows the belt's true
    ! per-axis origin after the union, which this subroutine does not
    ! compute) -- still a hard stop, now correct for every axis instead of
    ! y-only-relative-to-a-hardcoded-0.
    if (C_degen == 0.0d0) then
    do i = 1, ntotft
        ! checkIsOnFault (meshgen.f90) assumes a PLANAR vertical fault: one
        ! y value per fault, tested via fltxyz(1,2,iFault). A fault whose
        ! fymin /= fymax is not representable by that test at all (it would
        ! silently test against fymin only and drop every node past it).
        if (abs(fymax(i) - fymin(i)) > tol) then
            call abortRun(ERR_GEOM_MULTIFAULT_Y_BAD, &
                'checkInputConsistency: fault '//trim(itoa(i))//' has fymin /= fymax -- only a planar, ' // &
                'vertical fault (single y-plane) is supported; non-planar multi-fault geometry is out of scope.')
        endif
    enddo

    ! Distinct fault y-planes: two faults at the same y are not "two faults",
    ! they are one fault double-counted (and createMasterNode would build two
    ! overlapping split-node pairs at the same physical location).
    do i = 1, ntotft-1
        do j = i+1, ntotft
            if (abs(fymin(i) - fymin(j)) < tol) then
                call abortRun(ERR_GEOM_MULTIFAULT_Y_BAD, &
                    'checkInputConsistency: faults '//trim(itoa(i))//' and '//trim(itoa(j))// &
                    ' are both at y = '//trim(dtoa(fymin(i)))//' -- two faults must occupy distinct y-planes.')
            endif
        enddo
    enddo
    endif
end subroutine checkInputConsistency

function itoa(n) result(s)
    implicit none
    integer (kind = 4), intent(in) :: n
    character (len = 12) :: s
    write(s,'(I0)') n
end function itoa

function dtoa(x) result(s)
    use globalvar, only: dp
    implicit none
    real (kind = dp), intent(in) :: x
    character (len = 24) :: s
    write(s,'(F0.3)') x
end function dtoa

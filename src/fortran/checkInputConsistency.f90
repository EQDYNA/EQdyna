subroutine checkInputConsistency

    use globalvar
    use errorCodes
    implicit none
    integer (kind = 4) :: i, j
    real (kind = dp) :: yOverDy, nearestInt
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
    ! i<j distinctness/xz-match loops do not execute for ntotft<2).
    ! These multi-fault guards are about checkIsOnFault's PLANAR,
    ! C_degen==0 branch (one y-plane per fault, tested by exact-ish equality
    ! against fltxyz(1,2,iFault)) -- the narrow two-vertical-parallel-faults
    ! scope this release targets. C_degen>3 (wedge degeneration, e.g.
    ! test.tpv36/test.tpv37) is a DIFFERENT, pre-existing, already
    ! single-fault-only dipping-fault mechanism with a legitimate
    ! fymin /= fymax (the fault spans a range of grid y by construction) --
    ! orthogonal to this work and explicitly out of scope, so these guards
    ! do not apply to it.
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
        ! The y-grid's uniform belt spans [-dis4uniF*dy, +dis4uniB*dy] in
        ! exact steps of dy (meshgen.f90's getLocalOneDimCoorArrAndSize,
        ! dimId==2) -- a fault y not an integer multiple of dy falls between
        ! node lines and meshes with ZERO fault nodes, silently (eqquasi
        ! eeac6f9's exact failure mode: "Fault nodes = 0" read as cosmetic).
        yOverDy = fymin(i) / dy
        nearestInt = dble(nint(yOverDy))
        if (abs(yOverDy - nearestInt) > 1.0d-6) then
            call abortRun(ERR_GEOM_MULTIFAULT_Y_BAD, &
                'checkInputConsistency: fault '//trim(itoa(i))//'''s y = '//trim(dtoa(fymin(i)))// &
                ' is not an integer multiple of dy = '//trim(dtoa(dy))// &
                ' -- it would fall between mesh node lines and mesh with zero fault nodes.')
        endif
        if (fymin(i) < -dble(dis4uniF)*dy - tol .or. fymin(i) > dble(dis4uniB)*dy + tol) then
            call abortRun(ERR_GEOM_MULTIFAULT_Y_BAD, &
                'checkInputConsistency: fault '//trim(itoa(i))//'''s y = '//trim(dtoa(fymin(i)))// &
                ' lies outside the uniform-y mesh belt [-dis4uniF*dy, +dis4uniB*dy] = ['// &
                trim(dtoa(-dble(dis4uniF)*dy))//', '//trim(dtoa(dble(dis4uniB)*dy))// &
                ']; widen par.nuni_y_minus/par.nuni_y_plus or move the fault.')
        endif
        ! The uniform x/z belt (meshgen.f90's getLocalOneDimCoorArrAndSize,
        ! dimId==1/3) is anchored on fault 1's box ALONE (fltxyz(:,1,1),
        ! fltxyz(:,3,1)) -- independent per-fault x/z extents are out of
        ! scope for this release (two PARALLEL faults, same strike extent).
        if (i > 1) then
            if (abs(fxmin(i)-fxmin(1))>tol .or. abs(fxmax(i)-fxmax(1))>tol .or. &
                abs(fzmin(i)-fzmin(1))>tol .or. abs(fzmax(i)-fzmax(1))>tol) then
                call abortRun(ERR_GEOM_MULTIFAULT_XZ_BAD, &
                    'checkInputConsistency: fault '//trim(itoa(i))//' has a different x/z extent than fault 1. ' // &
                    'The shared uniform x/z mesh belt is built from fault 1 box alone, so every fault must ' // &
                    'share it (two parallel faults, same strike extent) -- independent per-fault x/z extents are out of scope.')
            endif
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

! Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
! MIT

! Subroutine list:
! 1. output_onfault_st
! 2. output_offfault_st
! 3. output_frt
! 4. output_timeanalysis
! 5. output_plastic_strain
! 6. find_surfNodeIdArr
! 7. output_gm
! 8. output_finalSurfDisp
! 9. output_src_evol
! 10. output_profile (+ readCpusAllowed, computeNumaNodes, parseRangeList,
!     readHostname, readPid, jsonNum helpers) -- docs/run_profile.md
! 11. report_dropped_offfault_st

!#1
subroutine output_onfault_st

    use globalvar
    implicit none
    
    integer (kind = 4) :: i, j 
    
    if(numOfOnFaultStCount>0) then
        do i=1,numOfOnFaultStCount
            j=anonfs(3,i)
            ! FIX (pathway_forward.md item 9): unit 51 used to be opened only
            ! when j==1 ("main fault stations") but written to unconditionally
            ! below -- for a station on fault 2+ (j>1) that skipped this whole
            ! block, so the write(51,...) a few lines down hit an unopened (or
            ! a previous station's already-closed) unit. This block is
            ! self-contained per loop iteration i already (open here, write
            ! below, close(51) at the bottom of this same iteration), so it
            ! must run for every station regardless of which fault j it is on
            ! -- there is nothing "main-fault-only" about building this
            ! station's own filename/header. No-op for the only case that
            ! currently runs (ntotft==1): there j is always 1, so this block
            ! always executed before too.
            sttmp = '    '
            dptmp = '   '
            ! Sign-aware strike field (item 22): (i3.3) overflowed to '***'
            ! for negative-x stations, colliding half of e.g. TPV29's list.
            ! SCEC convention: signed strike distance, zero-padded magnitude.
            if (nint(xonfs(1,anonfs(2,i),j)/100.d0) < 0) then
                write(sttmp,'(a1,i3.3)') '-', abs(nint(xonfs(1,anonfs(2,i),j)/100.d0))
            else
                write(sttmp,'(i3.3)') nint(xonfs(1,anonfs(2,i),j)/100.d0)
            endif
            write(dptmp,'(i3.3)') nint(abs(xonfs(2,anonfs(2,i),j))/dsin(fltxyz(2,4,1))/100.d0)
            open(51,file='faultst'//trim(adjustl(sttmp))//'dp'//trim(adjustl(dptmp))//'.txt',status='unknown')

            sttmp = '      '
            dptmp = '      '
            write(sttmp,'(f5.1)') xonfs(1,anonfs(2,i),j)/1000.d0
            ! Down-dip distance, not vertical depth: xonfs(2,...) is vertical
            ! depth (m, matched against nodeCoor(3) in meshgen.f90), so this
            ! must divide by sin(dip) to become a down-dip distance, exactly
            ! as the filename's own dp field does a few lines above (:48).
            ! Before this fix the stamp's value was unvalidated (never
            ! written) and equalled the filename's only because every
            ! currently-tested case with on-fault stations happens to have
            ! either dip=90 (sin=1) or C_degen==0 (fltxyz(2,4,*) forced to
            ! 90 regardless of the true mesh dip) -- test.tpv36/test.tpv37
            ! (C_degen=15) do have on-fault stations and would have printed
            ! vertical depth mislabeled "down-dip" had this shipped unfixed.
            write(dptmp,'(f5.1)') abs(xonfs(2,anonfs(2,i),j))/dsin(fltxyz(2,4,1))/1000.d0
            stLocStamp = '# location = on fault, '//trim(adjustl(sttmp))//' km along strike, '//trim(adjustl(dptmp))//' km down-dip'
            ! pathway item 85: stLocStamp was computed every call and never
            ! written -- pathway item 67's own evidence command looked for
            ! this exact line and found no file. Emit it as the header's
            ! first line.
            write(51,*) trim(stLocStamp)
            write(51,*) '# Project=',projectname
            write(51,*) '# Author=',author
            call date_and_time(values=dateTimeStamp)
            write(51,'( a10,i2,a1,i2,a1,i4,a1,i2,a1,i2,a1,i2)') ' # date = ',dateTimeStamp(2), &
                '/',dateTimeStamp(3),'/',dateTimeStamp(1),' ',dateTimeStamp(5),':',dateTimeStamp(6), &
                ':',dateTimeStamp(7)
            write(51,*) '# code = EQdyna'
            write(51,*) '# element_size =',dx
            write(51,'( a14,f8.4,a3)') '# time_step =', dt, ' s'
            write(51,'( a19,i6)') '# num_time_steps =', nstep
            write(51,*) '# Column #1 = Time (s)'
            write(51,*) '# Column #2 = horizontal slip (m)'
            write(51,*) '# Column #3 = horizontal slip rate (m/s)'
            write(51,*) '# Column #4 = horizontal shear stress (MPa)'
            write(51,*) '# Column #5 = down-dip slip (m)'
            write(51,*) '# Column #6 = down-dip slip rate (m/s)'
            write(51,*) '# Column #7 = down-dip shear stress (MPa)'
            write(51,*) '# Column #8 = normal stress (MPa)'
            if (friclaw>=3) then
                write(51,*) '# Column #9 = state variable psi (dimensionless)'
                write(51,*) '# Column #10 = Temperature (degrees Kelvin)'
                write(51,*) '# Column #11 = Pore pressure (MPa)'
                ! pathway item 67: this declaration used to be a single line,
                ! ABOVE this if-branch, unconditionally claiming 11 columns
                ! in "format E15.7" -- true for neither branch (this one
                ! writes E21.13 for column 1 and E16.7 for the rest, never
                ! E15.7). It now follows the branch and states this branch's
                ! true column count and its true, non-uniform format.
                write(51,*) '# Time series in 11 columns; column 1 in format E21.13, columns 2-11 in format E16.7'
                write(51,*) '# The line below lists the names of the data fields:'
                write(51,'(1X,103A)') 't h-slip h-slip-rate h-shear-stress v-slip v-slip-rate v-shear-stress n-stress psi temperature pressure'
                do j = 1, nstep
                    write(51,'( E21.13,10E16.7)') &
                        onFaultQuantHistSCECForm(1,j,i), &
                        onFaultQuantHistSCECForm(5,j,i), &
                        onFaultQuantHistSCECForm(2,j,i), &
                        onFaultQuantHistSCECForm(8,j,i)/1.0d6, &
                        -onFaultQuantHistSCECForm(6,j,i), &
                        -onFaultQuantHistSCECForm(3,j,i), &
                        -onFaultQuantHistSCECForm(9,j,i)/1.0d6, &
                        -onFaultQuantHistSCECForm(10,j,i)/1.0d6, &
                        onFaultQuantHistSCECForm(4,j,i), &
                        onFaultQuantHistSCECForm(12,j,i), &
                        onFaultQuantHistSCECForm(11,j,i)/1.0d6  
                enddo
            else
                ! pathway item 67: 8 columns on this branch (see the
                ! friclaw>=3 branch above for why this line moved here and
                ! why "E15.7" is gone -- the write below is E21.13 then
                ! 7x E16.7, matched here rather than the old unconditional,
                ! wrong-for-this-branch "11 columns in format E15.7").
                write(51,*) '# Time series in 8 columns; column 1 in format E21.13, columns 2-8 in format E16.7'
                write(51,*) '# The line below lists the names of the data fields:'
                write(51,'(1X,103A)') 't h-slip h-slip-rate h-shear-stress v-slip v-slip-rate v-shear-stress n-stress'
                do j = 1, nstep
                    write(51,'( E21.13,7E16.7)') &
                        onFaultQuantHistSCECForm(1,j,i), &
                        onFaultQuantHistSCECForm(5,j,i), &
                        onFaultQuantHistSCECForm(2,j,i), &
                        onFaultQuantHistSCECForm(8,j,i)/1.0d6, &
                        -onFaultQuantHistSCECForm(6,j,i), &
                        -onFaultQuantHistSCECForm(3,j,i), &
                        -onFaultQuantHistSCECForm(9,j,i)/1.0d6, &
                        -onFaultQuantHistSCECForm(10,j,i)/1.0d6
                enddo
            endif 
            close(51)
        enddo
    endif 
end subroutine output_onfault_st
!#2
subroutine output_offfault_st
    
    use globalvar
    implicit none
    
    integer (kind = 4) :: i, j     
    
    if(numOfOffFaultStCount>0) then
        do i=1,numOfOffFaultStCount
            bodytmp = '      '
            sttmp = '      '
            dptmp = '      '
            ! nint(), not int() (pathway item 93): the on-fault writer rounds
            ! and this one truncated, so a station at y=1.99 km filed as
            ! body010 instead of body020. The i4.3 edit carries the sign.
            write(bodytmp,'(i4.3)') nint(x4nds(2,OffFaultStNodeIdIndex(1,i))/100.d0)
            write(sttmp,'(i4.3)') nint(x4nds(1,OffFaultStNodeIdIndex(1,i))/100.d0)
            write(dptmp,'(i4.3)') nint(abs(x4nds(3,OffFaultStNodeIdIndex(1,i)))/100.d0)

            open(51,file='body'//trim(adjustl(bodytmp))//'st'//trim(adjustl(sttmp))//'dp'//trim(adjustl(dptmp))//'.txt',status='unknown')

            bodytmp = '      '
            sttmp = '      '
            dptmp = '      '
            write(bodytmp,'(f5.1)') x4nds(2,OffFaultStNodeIdIndex(1,i))/1000. 
            write(sttmp,'(f5.1)') x4nds(1,OffFaultStNodeIdIndex(1,i))/1000. 
            write(dptmp,'(f5.1)') abs(x4nds(3,OffFaultStNodeIdIndex(1,i)))/1000. 
            ! Vertical depth, NOT a down-dip distance -- unlike the on-fault
            ! stamp at :65, this one takes no /dsin(dip) correction. x4nds is
            ! read as a plain (x,y,z) triple in km from bStations.txt
            ! (readInputFiles.f90:218, scaled to m at :223) and its third
            ! component is matched against the raw mesh node z-coordinate
            ! nodeCoor(3) in setSurfaceStation (meshgen.f90:625). These are
            ! body/surface receivers at arbitrary (x,y,z) in the volume, off
            ! the fault plane by x4nds(2) -- a point that is not on the fault
            ! has no down-dip coordinate to convert to, so |z| is the right
            ! quantity and 'depth' is the right label.
            ! Missing ', ' separator (pathway item 88): the stamp used to
            ! render as '5.0 km along strike3.0 km depth'.
            stLocStamp = '# location = '//trim(adjustl(bodytmp))//' km off fault, '//trim(adjustl(sttmp))//' km along strike, '//trim(adjustl(dptmp))//' km depth'
            ! pathway item 88, same defect as item 85 one subroutine up:
            ! stLocStamp was computed every call and never written, so its
            ! value -- including the missing separator above -- was
            ! unverified. Emit it as the header's first line.
            write(51,*) trim(stLocStamp)
            write(51,*) '# Project=',projectname
            write(51,*) '# Author=',author
            call date_and_time(values=dateTimeStamp)
            write(51,'( a10,i2,a1,i2,a1,i4,a1,i2,a1,i2,a1,i2)') ' # date = ',dateTimeStamp(2), &
                    '/',dateTimeStamp(3),'/',dateTimeStamp(1),' ',dateTimeStamp(5),':',dateTimeStamp(6), &
                    ':',dateTimeStamp(7)
            write(51,*) '# code = EQdyna'
            write(51,*) '# element_size =',dx
            write(51,'( a14,f8.4,a3)') '# time_step=', dt, ' s'
            write(51,'( a19,i6)') '# num_time_steps=', nstep
            ! pathway item 88 (item 67 class): this legend used to carry a
            ! spurious second '# Column #3 = horizontal displacement (m)'
            ! line, so it declared 8 entries for a 7-column write and every
            ! name from #3 down sat against the wrong ordinal. The list below
            ! is derived from the write(51,'( E21.13,6E16.7)') statement at
            ! the bottom of this loop, via the idhist row order built in
            ! eqdyna3d.f90:180-190 (per station: iDof=1..3 outer,
            ! disp/vel inner), which fixes OffFaultStGramSCEC row
            ! (i-1)*6+1+2*(iDof-1)+dispOrVel:
            !   +2 = x disp, +3 = x vel, +4 = y disp,
            !   +5 = y vel,  +6 = z disp, +7 = z vel.
            ! Write order is therefore t, x-disp, x-vel, -z-disp, -z-vel,
            ! y-disp, y-vel. x is along strike ('horizontal'), y is
            ! fault-normal ('normal'), and z is negated because the model's
            ! z axis points up while the SCEC 'vertical' component is
            ! positive down -- the same negation the on-fault writer applies
            ! to its down-dip components at :109-112.
            ! Item 93: declare the column count and format like the on-fault
            ! writer does, derived from the write(51,'( E21.13,6E16.7)')
            ! below -- item 67's class in its absent form.
            write(51,*) '# Time series in 7 columns; column 1 in format E21.13, columns 2-7 in format E16.7'
            write(51,*) '# Column #1 = Time (s)'
            write(51,*) '# Column #2 = horizontal displacement (m)'
            write(51,*) '# Column #3 = horizontal velocity (m/s)'
            write(51,*) '# Column #4 = vertical displacement (m)'
            write(51,*) '# Column #5 = vertical velocity (m/s)'
            write(51,*) '# Column #6 = normal displacement (m)'
            write(51,*) '# Column #7 = normal velocity (m/s)'
            write(51,*) '#'
            write(51,*) '# The line below lists the names of the data fields:'
            ! Item 93: the same '(1X,103A)' edit as output_onfault_st, not
            ! list-directed, so both station writers emit one header shape.
            write(51,'(1X,103A)') 't h-disp h-vel v-disp v-vel n-disp n-vel'
            do j=1, nstep 
                write(51,'( E21.13,6E16.7)') &
                    OffFaultStGramSCEC(1,j), &
                    OffFaultStGramSCEC((i-1)*6+2,j), &
                    OffFaultStGramSCEC((i-1)*6+3,j), &
                    -OffFaultStGramSCEC((i-1)*6+6,j), &
                    -OffFaultStGramSCEC((i-1)*6+7,j), &
                    OffFaultStGramSCEC((i-1)*6+4,j), &
                    OffFaultStGramSCEC((i-1)*6+5,j)
            enddo
            close(51)
        enddo
    endif
end subroutine output_offfault_st

!#3
subroutine output_frt
    ! The subroutine output_frt generates frt.txt* files for each MPI process.
    ! frt.txt* contain on-fault variables for visualization 
    !   and restart files for the next deformation phase. 
    use globalvar
    implicit none
    
    integer (kind = 4) :: i, j
    integer (kind = 4), parameter :: UNIT_FRT_BASE = 10004

    if(nftnd(1) > 0) then
        open(unit=UNIT_FRT_BASE+me,file='frt.txt'//mm,status='unknown')

        write(UNIT_FRT_BASE+me,'(1x,22e18.7e4)')    &
                ! 3 coordinates of the fault nodes.
            ((meshCoor(j,nsmp(1,i,1)), j = 1,3), &
                ! rupture time
            fnft(i,1),                    &
                ! 71-73: final slips, slipd, slipn
                ! 74-76: final sliprates, sliprated, slipraten
            (fric(j,i,1), j = FRIC_SLOT_SLIP_STRIKE, FRIC_SLOT_SLIPRATE_MAX),     &
                ! 47: final slip rate
            fric(FRIC_SLOT_PEAK_SLIPRATE,i,1),                 &
                ! 78: final effective normal stress, tnrm
                ! 79: final shear strike stress, tstk
                ! 80: final shear dip stress, tdip
            fric(FRIC_SLOT_TRACT_NORM,i,1),                 &
            fric(FRIC_SLOT_TRACT_STRIKE,i,1),                 &
            fric(FRIC_SLOT_TRACT_DIP,i,1),                 &
                ! 31-33, vxm, vym, vzm, 3 vel components of master nodes.
                ! 34-36, vxs, vys, vzs, 3 vel components of slave nodes.
            fric(FRIC_SLOT_VEL_MASTER_X,i,1),                 &
            fric(FRIC_SLOT_VEL_MASTER_Y,i,1),                 &
            fric(FRIC_SLOT_VEL_MASTER_Z,i,1),                 &
            fric(FRIC_SLOT_VEL_SLAVE_X,i,1),                 &
            fric(FRIC_SLOT_VEL_SLAVE_Y,i,1),                 &
            fric(FRIC_SLOT_VEL_SLAVE_Z,i,1),                 &
                ! 20: state variable in RSF
            fric(FRIC_SLOT_STATE,i,1),                 &
                ! 21: state variable for normal stress variation (Shi and Day)
            fric(FRIC_SLOT_THETA_PC,i,1),                 &
                !
            i=1,nftnd(1)) ! Finish the write(10004,me, ...) line.

        close(UNIT_FRT_BASE+me)
    endif
end subroutine output_frt

!#4
subroutine output_timeanalysis

    use globalvar
    implicit none
    
    integer (kind = 4) :: i, j     
    
    open(unit=14,file='compTime'//mm,status='unknown')    !rupture time
        write(14,'(1x,10e18.7e4,2i10)') (compTimeInSeconds(i),i=1,9),MPICommTimeInSeconds,totalNumOfElements,totalNumOfEquations
    close(14)
end subroutine output_timeanalysis

!#5
subroutine output_plastic_strain
    use globalvar
    implicit none
    integer (kind = 4) :: i, j
    integer (kind = 4), parameter :: UNIT_PSTR_BASE = 10007
    real (kind = dp) :: sc(3)
    if (output_plastic == 1) then

        ! The output window is par.plasticOutputHalfWidth, read from
        ! bGlobal.txt. It defaults to (5, 2, 8) km, the box this subroutine
        ! hardcoded until v5.9.0 -- which was sized for one case and would
        ! silently clip, or entirely miss, a case with a different fault or
        ! domain extent (TPV30's fault alone is 40 x 20 km).
        do i=1,totalNumOfElements
                    if ((pstrain(i)>1.0d-4).and.(abs(meshCoor(1,nodeElemIdRelation(1,i)))<plasticOutputHalfWidth(1)).and.(abs(meshCoor(2,nodeElemIdRelation(1,i)))<plasticOutputHalfWidth(2)).and.(abs(meshCoor(3,nodeElemIdRelation(1,i)))<plasticOutputHalfWidth(3))) then
                open(unit=UNIT_PSTR_BASE+me,file='pstr.txt'//mm,status='unknown',position='append')
                sc=0.0d0
                do j=1,8
                    sc(1)=sc(1)+meshCoor(1,nodeElemIdRelation(j,i))
                    sc(2)=sc(2)+meshCoor(2,nodeElemIdRelation(j,i))
                    sc(3)=sc(3)+meshCoor(3,nodeElemIdRelation(j,i))
                enddo
                sc(1)=sc(1)/8.0d0
                sc(2)=sc(2)/8.0d0
                sc(3)=sc(3)/8.0d0            
                write(UNIT_PSTR_BASE+me,'(1x,16e18.7e4)') sc(1),sc(2),sc(3),pstrain(i),(stressArr(stressCompIndexArr(i)+j),j=1,12)
            endif
        enddo
    endif
end subroutine output_plastic_strain

!#6
subroutine find_surfaceNodeIdArr
    use globalvar
    implicit none
    integer (kind = 4) :: i, j
    integer (kind = 4), parameter :: UNIT_SURFCOOR_BASE = 10008
    real (kind = dp), parameter :: STATION_SEARCH_HALFWIDTH_M = 20.0d3 ! along-strike/along-strike-normal search box half-width around the fault trace, m
    real (kind = dp) :: sc(3)
    if (outputGroundMotion==1 .or. outputFinalSurfDisp==1) then
        do i=1,totalNumOfNodes
            if ((meshCoor(1,i)<fltxyz(2,1,1)+STATION_SEARCH_HALFWIDTH_M) .and. (meshCoor(1,i)>fltxyz(1,1,1)-STATION_SEARCH_HALFWIDTH_M) &
                    .and. (meshCoor(2,i)<fltxyz(2,2,1)+STATION_SEARCH_HALFWIDTH_M) .and. (meshCoor(2,i)>fltxyz(1,2,1)-STATION_SEARCH_HALFWIDTH_M) &
                    .and. (abs(meshCoor(3,i))<dx/1000)) then
                surface_nnode = surface_nnode + 1
                surfaceNodeIdArr(surface_nnode) = i
                open(unit=UNIT_SURFCOOR_BASE+me,file='surface_coor.txt'//mm,status='unknown',position='append')
                    write(UNIT_SURFCOOR_BASE+me,'(1x,3e18.7e4)') meshCoor(1,i), meshCoor(2,i), meshCoor(3,i)
                close(UNIT_SURFCOOR_BASE+me)
            endif
        enddo
    endif
end subroutine find_surfaceNodeIdArr

!#7
subroutine output_gm
    use globalvar
    implicit none
    integer (kind = 4) :: i, j, nodeId
    integer (kind = 4), parameter :: UNIT_GM_BASE = 10009

    if (outputGroundMotion == 1 .and. surface_nnode > 0) then
        open(unit=UNIT_GM_BASE+me,file='gm'//mm,status='unknown',position='append', access='stream')
            do i=1,surface_nnode
                nodeId = surfaceNodeIdArr(i)
                write(UNIT_GM_BASE+me) velArr(1,nodeId), velArr(2,nodeId), velArr(3,nodeId)
            enddo
    endif
end subroutine output_gm

!#8
subroutine output_finalSurfDisp
    use globalvar
    implicit none
    integer (kind = 4) :: i, j, nodeId
    integer (kind = 4), parameter :: UNIT_FINALSURFDISP_BASE = 20009

    if (outputFinalSurfDisp == 1 .and. surface_nnode > 0) then
        open(unit=UNIT_FINALSURFDISP_BASE+me,file='finalSurfDisp.txt'//mm, status='unknown')
            do i=1,surface_nnode
                nodeId = surfaceNodeIdArr(i)
                write(UNIT_FINALSURFDISP_BASE+me, '(1x,3e18.7e4)') dispArr(1,nodeId), dispArr(2,nodeId), dispArr(3,nodeId)
            enddo
    endif
end subroutine output_finalSurfDisp

!#9
subroutine output_src_evol
    ! The subroutine output_src_evol generates binary src_evol files for each MPI process.
    ! srv_evol contains on-fault slip-rate system states for AI and visualization. 
    use globalvar
    implicit none
    
    integer (kind = 4) :: i, j, nodeId
    integer (kind = 4), parameter :: UNIT_SRC_EVOL_BASE = 30009

    if(nftnd(1) > 0) then
        open(unit=UNIT_SRC_EVOL_BASE+me,file='src_evol'//mm,position='append', access='stream')
            do i=1,nftnd(1)
                write(UNIT_SRC_EVOL_BASE+me) fric(FRIC_SLOT_PEAK_SLIPRATE,i,1) ! 47: final slip rate
            enddo
        close(UNIT_SRC_EVOL_BASE+me)
    endif
end subroutine output_src_evol

!#10
subroutine output_profile(setupBucket, elementBucket, faultBucket, exchangeBucket, &
                           waitBucket, ioBucket, loopS, totalS, nstepsArg, samplingEvery)
    ! ALWAYS-ON per-rank profile.rank<r>.json (docs/run_profile.md). Off
    ! switch EQDYNA_PROFILE=0 exists only for the overhead A/B; default on.
    ! Every value here is a directly-measured quantity or a difference of
    ! two directly-measured monotonic counters at named checkpoints --
    ! unaccounted_s is a REPORTED diagnostic, never fed back into a bucket.
    use globalvar
    implicit none
    real (kind = dp), intent(in) :: setupBucket, elementBucket, faultBucket, exchangeBucket, &
                                     waitBucket, ioBucket, loopS, totalS
    integer (kind = 4), intent(in) :: nstepsArg, samplingEvery
    integer (kind = 4), parameter :: UNIT_PROFILE_BASE = 40009
    integer (kind = 4), allocatable :: cpuList(:), numaList(:)
    integer (kind = 4) :: nCpus, nNuma, i, pid
    character(len=256) :: hostStr
    real (kind = dp) :: unaccounted

    ! EQDYNA_PROFILE is read ONCE, at startup (eqdyna3d.f90, right after
    ! MPI_Init), into the module logical profileEnabled (globalvar.f90).
    ! Re-reading it here via get_environment_variable broke that read-once
    ! contract AND used an 8-char buffer that silently truncated anything
    ! longer -- both defects are avoided by testing the already-parsed
    ! module flag instead.
    if (.not. profileEnabled) return

    call readCpusAllowed(cpuList, nCpus)
    call computeNumaNodes(cpuList, nCpus, numaList, nNuma)
    call readHostname(hostStr)
    call readPid(pid)
    unaccounted = totalS - (setupBucket + elementBucket + faultBucket + exchangeBucket + waitBucket + ioBucket)

    open(unit=UNIT_PROFILE_BASE+me, file='profile.rank'//trim(mm)//'.json', status='unknown')
    write(UNIT_PROFILE_BASE+me,'(A)') '{'
    write(UNIT_PROFILE_BASE+me,'(A)') '"schema": "eqdyna-profile/1",'
    write(UNIT_PROFILE_BASE+me,'(A)') '"backend": "fortran",'
    write(UNIT_PROFILE_BASE+me,'(A,I0,A)') '"rank": ', me, ','
    write(UNIT_PROFILE_BASE+me,'(A,I0,A)') '"nranks": ', totalNumOfMPIProcs, ','
    write(UNIT_PROFILE_BASE+me,'(A,A,A)') '"host": "', trim(hostStr), '",'
    write(UNIT_PROFILE_BASE+me,'(A,I0,A)') '"pid": ', pid, ','
    write(UNIT_PROFILE_BASE+me,'(A)',advance='no') '"cpus_allowed": ['
    do i = 1, nCpus
        if (i > 1) write(UNIT_PROFILE_BASE+me,'(A)',advance='no') ','
        write(UNIT_PROFILE_BASE+me,'(I0)',advance='no') cpuList(i)
    enddo
    write(UNIT_PROFILE_BASE+me,'(A)') '],'
    write(UNIT_PROFILE_BASE+me,'(A)',advance='no') '"numa_nodes": ['
    do i = 1, nNuma
        if (i > 1) write(UNIT_PROFILE_BASE+me,'(A)',advance='no') ','
        write(UNIT_PROFILE_BASE+me,'(I0)',advance='no') numaList(i)
    enddo
    write(UNIT_PROFILE_BASE+me,'(A)') '],'
    write(UNIT_PROFILE_BASE+me,'(A,I0,A)') '"nsteps": ', nstepsArg, ','
    write(UNIT_PROFILE_BASE+me,'(A,I0,A)') '"sampling_every": ', samplingEvery, ','
    write(UNIT_PROFILE_BASE+me,'(A)') '"buckets_s": {'
    write(UNIT_PROFILE_BASE+me,'(A,A,A)') '"setup": ', trim(jsonNum(setupBucket)), ','
    write(UNIT_PROFILE_BASE+me,'(A,A,A)') '"element": ', trim(jsonNum(elementBucket)), ','
    write(UNIT_PROFILE_BASE+me,'(A,A,A)') '"fault": ', trim(jsonNum(faultBucket)), ','
    write(UNIT_PROFILE_BASE+me,'(A,A,A)') '"exchange": ', trim(jsonNum(exchangeBucket)), ','
    write(UNIT_PROFILE_BASE+me,'(A,A,A)') '"wait": ', trim(jsonNum(waitBucket)), ','
    write(UNIT_PROFILE_BASE+me,'(A,A)') '"io": ', trim(jsonNum(ioBucket))
    write(UNIT_PROFILE_BASE+me,'(A)') '},'
    write(UNIT_PROFILE_BASE+me,'(A,A,A)') '"loop_s": ', trim(jsonNum(loopS)), ','
    write(UNIT_PROFILE_BASE+me,'(A,A,A)') '"total_s": ', trim(jsonNum(totalS)), ','
    write(UNIT_PROFILE_BASE+me,'(A,A)') '"unaccounted_s": ', trim(jsonNum(unaccounted))
    write(UNIT_PROFILE_BASE+me,'(A)') '}'
    close(UNIT_PROFILE_BASE+me)
contains
    function jsonNum(x) result(s)
        real (kind = dp), intent(in) :: x
        character(len=32) :: s
        write(s,'(F24.9)') x
        s = adjustl(s)
    end function jsonNum

    subroutine readCpusAllowed(cpuListOut, nCpusOut)
        ! This process's CPU affinity, from /proc/self/status
        ! "Cpus_allowed_list:" (constraint 3: read via Fortran, not faked,
        ! not shelled out). Internal procedure: allocatable dummies need an
        ! explicit interface, which an internal procedure gets for free.
        integer, allocatable, intent(out) :: cpuListOut(:)
        integer, intent(out) :: nCpusOut
        integer :: u, ios, colonPos
        character(len=512) :: line
        logical :: found

        found = .false.
        open(newunit=u, file='/proc/self/status', status='old', action='read', iostat=ios)
        if (ios /= 0) then
            allocate(cpuListOut(0)); nCpusOut = 0; return
        endif
        do
            read(u,'(A)',iostat=ios) line
            if (ios /= 0) exit
            if (index(line,'Cpus_allowed_list:') == 1) then
                colonPos = index(line, ':')
                call parseRangeList(line(colonPos+1:), cpuListOut, nCpusOut)
                found = .true.
                exit
            endif
        enddo
        close(u)
        if (.not. found) then
            allocate(cpuListOut(0)); nCpusOut = 0
        endif
    end subroutine readCpusAllowed

    subroutine computeNumaNodes(cpuListIn, nCpusIn, numaListOut, nNumaOut)
        ! NUMA node(s) backing cpuListIn, from
        ! /sys/devices/system/node/node<k>/cpulist (constraint 3). Probes
        ! node indices 0..31 by open() success/failure -- there is no
        ! portable Fortran directory listing.
        integer, intent(in) :: cpuListIn(:), nCpusIn
        integer, allocatable, intent(out) :: numaListOut(:)
        integer, intent(out) :: nNumaOut
        integer :: k, u, ios, i, j
        character(len=256) :: fname
        character(len=512) :: line
        integer, allocatable :: nodeCpus(:), tmp(:)
        integer :: nNodeCpus
        logical :: overlap

        allocate(tmp(64))
        nNumaOut = 0
        do k = 0, 31
            write(fname,'(A,I0,A)') '/sys/devices/system/node/node', k, '/cpulist'
            open(newunit=u, file=trim(fname), status='old', action='read', iostat=ios)
            if (ios /= 0) cycle
            read(u,'(A)',iostat=ios) line
            close(u)
            if (ios /= 0) cycle
            call parseRangeList(line, nodeCpus, nNodeCpus)
            overlap = .false.
            do i = 1, nCpusIn
                do j = 1, nNodeCpus
                    if (cpuListIn(i) == nodeCpus(j)) then
                        overlap = .true.
                        exit
                    endif
                enddo
                if (overlap) exit
            enddo
            if (overlap) then
                nNumaOut = nNumaOut + 1
                tmp(nNumaOut) = k
            endif
            if (allocated(nodeCpus)) deallocate(nodeCpus)
        enddo
        allocate(numaListOut(nNumaOut))
        numaListOut(1:nNumaOut) = tmp(1:nNumaOut)
        deallocate(tmp)
    end subroutine computeNumaNodes

    subroutine parseRangeList(line, arrOut, nOut)
        ! Shared parser for kernel-style cpu range lists, e.g.
        ! "48-51,60,62-63". Used for both /proc/self/status
        ! Cpus_allowed_list and /sys/devices/system/node/node<k>/cpulist,
        ! which share this format.
        character(len=*), intent(in) :: line
        integer, allocatable, intent(out) :: arrOut(:)
        integer, intent(out) :: nOut
        character(len=len(line)) :: buf
        integer :: i, n, lo, hi, dashPos, pos, startTok
        integer, allocatable :: tmp(:)

        buf = trim(adjustl(line))
        n = len_trim(buf)
        allocate(tmp(4096))
        nOut = 0
        startTok = 1
        i = 1
        do while (i <= n+1)
            if (i > n .or. buf(i:i) == ',') then
                if (i > startTok) then
                    dashPos = index(buf(startTok:i-1), '-')
                    if (dashPos > 0) then
                        read(buf(startTok:startTok+dashPos-2), *) lo
                        read(buf(startTok+dashPos:i-1), *) hi
                    else
                        read(buf(startTok:i-1), *) lo
                        hi = lo
                    endif
                    do pos = lo, hi
                        nOut = nOut + 1
                        tmp(nOut) = pos
                    enddo
                endif
                startTok = i + 1
            endif
            i = i + 1
        enddo
        allocate(arrOut(nOut))
        arrOut(1:nOut) = tmp(1:nOut)
        deallocate(tmp)
    end subroutine parseRangeList
end subroutine output_profile

subroutine readHostname(hostStr)
    ! /proc/sys/kernel/hostname, read directly (no shell-out).
    implicit none
    character(len=256), intent(out) :: hostStr
    integer :: u, ios
    hostStr = 'unknown'
    open(newunit=u, file='/proc/sys/kernel/hostname', status='old', action='read', iostat=ios)
    if (ios == 0) then
        read(u,'(A)',iostat=ios) hostStr
        close(u)
    endif
end subroutine readHostname

subroutine readPid(pid)
    ! /proc/self/status "Pid:" line, same file and technique as
    ! output_profile's internal readCpusAllowed, per constraint 3 (Fortran affinity via /proc).
    implicit none
    integer, intent(out) :: pid
    integer :: u, ios, colonPos
    character(len=256) :: line
    pid = -1
    open(newunit=u, file='/proc/self/status', status='old', action='read', iostat=ios)
    if (ios /= 0) return
    do
        read(u,'(A)',iostat=ios) line
        if (ios /= 0) exit
        if (index(line,'Pid:') == 1) then
            colonPos = index(line, ':')
            read(line(colonPos+1:), *) pid
            exit
        endif
    enddo
    close(u)
end subroutine readPid

!#11
subroutine report_dropped_offfault_st(matchedAnyRank)
    ! pathway_forward.md item 94. setSurfaceStation (meshgen.f90) matches an
    ! off-fault station's DEPTH exactly (|nodeCoor(3) - x4nds(3,i)| < tol)
    ! and snaps only x and y to the nearest interior node, so a station whose
    ! depth is not a grid z-plane, or which lies outside the interior x/y
    ! range, matches no node on any rank and output_offfault_st writes no
    ! file for it. That used to happen with no message at all (test.tpv8 at
    ! dx = 500 m: its four z = -0.3 km stations). This names every such
    ! station. It does NOT snap or refuse: either would change which body*
    ! files a case writes, and that is the owner's call.
    ! matchedAnyRank(i) is (n4yn(i) /= 0) OR-reduced over all ranks
    ! (checkOffFaultStationCoverage, eqdyna3d.f90), so a station found by any
    ! rank counts as written.
    use globalvar
    implicit none

    logical, intent(in) :: matchedAnyRank(totalNumOfOffSt)
    integer (kind = 4) :: i, nDropped

    nDropped = count(.not. matchedAnyRank)
    if (nDropped == 0) return

    write(*,'(a,i0,a,i0,a)') ' WARNING: ', nDropped, ' of ', totalNumOfOffSt, &
        ' requested off-fault stations match no grid node and get NO body* file'
    write(*,'(a)') '   (setSurfaceStation, meshgen.f90: depth must equal a grid z-plane' // &
        ' within tol; x and y snap to the nearest interior node)'
    do i = 1, totalNumOfOffSt
        if (.not. matchedAnyRank(i)) then
            write(*,'(a,i0,a,3f10.3,a)') '   dropped off-fault station ', i, &
                ' at x,y,z =', x4nds(1,i)/1000.d0, x4nds(2,i)/1000.d0, x4nds(3,i)/1000.d0, ' km'
        endif
    enddo
end subroutine report_dropped_offfault_st

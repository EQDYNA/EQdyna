#!/usr/bin/env python3
"""Apply the row-149 force-budget instrumentation on top of the earlier diagTrace patch.
usage: apply_budget_edits.py <tree_root>"""
import sys, pathlib
root = pathlib.Path(sys.argv[1]) / 'src' / 'fortran'

def rep(fn, old, new, count=1):
    p = root / fn
    s = p.read_text()
    n = s.count(old)
    if n != count:
        raise SystemExit(f'{fn}: expected {count} occurrence(s), found {n}: {old[:60]!r}')
    s2 = s.replace(old, new)
    tmp = p.with_suffix('.tmp')
    tmp.write_text(s2)
    tmp.replace(p)

# ---- globalvar.f90
rep('globalvar.f90',
"    real (kind = dp) :: rdampm=0.0d0, rdampk        ! Rayleigh mass-/stiffness-proportional damping coefficients\n",
"    real (kind = dp) :: rdampm=0.0d0, rdampk        ! Rayleigh mass-/stiffness-proportional damping coefficients\n"
"    ! row 149 force-budget diagnostic (scratch build only, never shipped)\n"
"    integer (kind = 4) :: nDiagEl = 0, nDiagPr = 0\n"
"    integer (kind = 4), allocatable :: diagElIdx(:), diagNodeFlag(:), diagPrList(:)\n"
"    integer (kind = 8), allocatable :: diagElList(:)\n"
"    real (kind = dp), allocatable :: diagElF(:,:,:,:), diagElS(:,:), diagFpre(:)\n"
"    real (kind = dp) :: lastRjust = 1.0d0\n")

# ---- calcElemKU.f90
rep('calcElemKU.f90', "    pstrmag = 0.0d0\n    stressrate = 0.0d0\n",
    "    pstrmag = 0.0d0\n    lastRjust = 1.0d0\n    stressrate = 0.0d0\n")
rep('calcElemKU.f90', "            rjust=yield/taomax + (1-yield/taomax)*exp(-dt/tv)\n",
    "            rjust=yield/taomax + (1-yield/taomax)*exp(-dt/tv)\n            lastRjust = rjust\n")

# ---- assembleGlobalKU.f90
rep('assembleGlobalKU.f90',
"    integer (kind = 4) :: i, j\n    real (kind = dp) :: pstrinc, efPML(96), elresf(nee), al(ned,nen)\n",
"    integer (kind = 4) :: i, j, idg\n    real (kind = dp) :: pstrinc, efPML(96), elresf(nee), al(ned,nen), fBtS(nee), gacc\n")
rep('assembleGlobalKU.f90',
"        call calcElemMass(elemass(1:nee,nel),al,elresf)\n",
"        call calcElemMass(elemass(1:nee,nel),al,elresf)\n"
"        idg = 0\n"
"        if (allocated(diagElIdx)) idg = diagElIdx(nel)\n"
"        if (idg > 0) then\n"
"            ! parts 1 (gravity body force) and 2 (mass-proportional damping)\n"
"            gacc = (1.0d0-C_elastic)*grav*(roumax-(gamar+1.0d0)*rhow)/roumax\n"
"            diagElF(:,:,:,idg) = 0.0d0\n"
"            do i = 1, nen\n"
"                do j = 1, ned\n"
"                    diagElF(j,i,2,idg) = -rdampm*velArr(j,nodeElemIdRelation(i,nel))*elemass((i-1)*ned+j,nel)\n"
"                enddo\n"
"                diagElF(3,i,1,idg) = -gacc*elemass((i-1)*ned+3,nel)\n"
"            enddo\n"
"        endif\n")
rep('assembleGlobalKU.f90',
"            pstrain(nel) = pstrain(nel) + pstrinc\n",
"            pstrain(nel) = pstrain(nel) + pstrinc\n"
"            if (idg > 0) then\n"
"                ! part 3: B^T sigma (stress divergence); part 4: stiffness-proportional damping (remainder)\n"
"                call diagBtS(eleshp(1,1,nel), stressArr(stressCompIndexArr(nel)+1:stressCompIndexArr(nel)+6), -eledet(nel), fBtS)\n"
"                do i = 1, nen\n"
"                    do j = 1, ned\n"
"                        diagElF(j,i,3,idg) = fBtS((i-1)*ned+j)\n"
"                        diagElF(j,i,4,idg) = elresf((i-1)*ned+j) - diagElF(j,i,1,idg) - diagElF(j,i,2,idg) - fBtS((i-1)*ned+j)\n"
"                    enddo\n"
"                enddo\n"
"                diagElS(1:6,idg) = stressArr(stressCompIndexArr(nel)+1:stressCompIndexArr(nel)+6)\n"
"                diagElS(7,idg) = lastRjust\n"
"                diagElS(8,idg) = pstrain(nel)\n"
"                diagElS(9,idg) = eledet(nel)*w\n"
"            endif\n")
rep('assembleGlobalKU.f90',
"end subroutine assembleGlobalKU\n",
"end subroutine assembleGlobalKU\n\n"
"subroutine diagBtS(globalShapeFunc, stress, constk, fout)\n"
"    ! row 149 diagnostic: the B^T*sigma nodal force of calcElemKU, stress part only\n"
"    use globalvar\n"
"    implicit none\n"
"    integer (kind = 4) :: i, j1, j2, j3\n"
"    real (kind = dp) :: globalShapeFunc(nrowsh-1,nen), stress(6), constk, fout(nee), bb(nrowb,nee), st(6)\n"
"    call calcB(globalShapeFunc,bb)\n"
"    st = constk*w*stress\n"
"    do i=1,nen\n"
"        j1 = ned * (i - 1) + 1\n"
"        j2 = ned * (i - 1) + 2\n"
"        j3 = ned * (i - 1) + 3\n"
"        fout(j1) = bb(1,j1)*st(1) + bb(5,j1)*st(5) + bb(6,j1)*st(6)\n"
"        fout(j2) = bb(2,j2)*st(2) + bb(4,j2)*st(4) + bb(6,j2)*st(6)\n"
"        fout(j3) = bb(3,j3)*st(3) + bb(4,j3)*st(4) + bb(5,j3)*st(5)\n"
"    enddo\n"
"end subroutine diagBtS\n")

# ---- calcHourglassResist.f90
rep('calcHourglassResist.f90',
"    integer (kind = 4) :: i , j, itmp, fi(4,8)\n",
"    integer (kind = 4) :: i , j, itmp, fi(4,8), idg\n    real (kind = dp) :: phide(ned), fe(ned)\n")
rep('calcHourglassResist.f90',
"    do nel = 1, totalNumOfElements\n        do i = 1, nen\n            do j = 1, ned\n                vl(j,i)",
"    do nel = 1, totalNumOfElements\n        idg = 0\n        if (allocated(diagElIdx)) idg = diagElIdx(nel)\n        do i = 1, nen\n            do j = 1, ned\n                vl(j,i)")
rep('calcHourglassResist.f90',
"                !... assemble to global right-hand-side force vector\n",
"                if (idg > 0) then\n"
"                    ! part 5: hourglass force on displacement; part 6: its rdampk*velocity part\n"
"                    do i = 1,ned\n"
"                        phide(i) = 0.0d0\n"
"                        do j = 1,nen\n"
"                            phide(i) = phide(i) + phi(j,itmp,nel) * dispArr(i,nodeElemIdRelation(j,nel))\n"
"                        enddo\n"
"                    enddo\n"
"                    do i = 1,nen\n"
"                        fe(1) = phi(i,itmp,nel) * (ss(1,nel)*phide(1) + ss(2,nel)*phide(2) + ss(3,nel)*phide(3))\n"
"                        fe(2) = phi(i,itmp,nel) * (ss(2,nel)*phide(1) + ss(4,nel)*phide(2) + ss(5,nel)*phide(3))\n"
"                        fe(3) = phi(i,itmp,nel) * (ss(3,nel)*phide(1) + ss(5,nel)*phide(2) + ss(6,nel)*phide(3))\n"
"                        diagElF(1:3,i,5,idg) = diagElF(1:3,i,5,idg) - fe(1:3)\n"
"                        diagElF(1:3,i,6,idg) = diagElF(1:3,i,6,idg) - (fhr(1:3,i) - fe(1:3))\n"
"                    enddo\n"
"                endif\n"
"                !... assemble to global right-hand-side force vector\n")
rep('calcHourglassResist.f90',
"            do i = 1, nen\n                do j = 1, ned\n                    if (numOfDofPerNodeArr(nodeElemIdRelation(i,nel)).eq.3) then\n",
"            if (idg > 0) then\n"
"                do i = 1, nen\n"
"                    do j = 1, ned\n"
"                        diagElF(j,i,6,idg) = f((i-1)*3+j)\n"
"                    enddo\n"
"                enddo\n"
"            endif\n"
"            do i = 1, nen\n                do j = 1, ned\n                    if (numOfDofPerNodeArr(nodeElemIdRelation(i,nel)).eq.3) then\n")

# ---- driver.f90 (on top of the earlier diagTrace patch)
rep('driver.f90', "        nodalForceArr = 0.0d0\n        \n        call assembleGlobalKU\n",
    "        nodalForceArr = 0.0d0\n        \n        if (nt == 1) call diagBudgetInit\n        call assembleGlobalKU\n")
rep('driver.f90', "        call diagTrace\n        if (friclaw == 5)", "        call diagTrace\n        call diagBudget(1)\n        if (friclaw == 5)")
rep('driver.f90', "        call faulting\n        nodalForceArr(1:totalNumOfEquations) = nodalForceArr(1:totalNumOfEquations)/nodalMassArr(1:totalNumOfEquations)\n",
    "        call faulting\n        call diagBudget(2)\n        nodalForceArr(1:totalNumOfEquations) = nodalForceArr(1:totalNumOfEquations)/nodalMassArr(1:totalNumOfEquations)\n")

INIT_SUB = r'''
subroutine diagBudgetInit
! row 149 force-budget diagnostic: select the split-node pairs at x = 0 in the
! top three fault rows and every element touching one of them.
    use globalvar
    implicit none
    integer (kind = 4) :: i, k, ift
    integer (kind = 8) :: nel, ns, nm
    character (len = 64) :: fname
    ift = 1
    allocate(diagElIdx(totalNumOfElements), diagNodeFlag(totalNumOfNodes), diagPrList(max(nftnd(ift),1)))
    allocate(diagFpre(totalNumOfEquations))
    diagElIdx = 0
    diagNodeFlag = 0
    diagPrList = 0
    diagFpre = 0.0d0
    nDiagPr = 0
    do i = 1, nftnd(ift)
        ns = nsmp(1,i,ift)
        nm = nsmp(2,i,ift)
        if (abs(meshCoor(1,ns)) > 1.0d0) cycle
        if (meshCoor(3,ns) < -2.5d0*dz) cycle
        nDiagPr = nDiagPr + 1
        diagPrList(nDiagPr) = i
        diagNodeFlag(ns) = 1
        diagNodeFlag(nm) = 2
    enddo
    nDiagEl = 0
    do nel = 1, totalNumOfElements
        do k = 1, nen
            if (diagNodeFlag(nodeElemIdRelation(k,nel)) > 0) then
                nDiagEl = nDiagEl + 1
                diagElIdx(nel) = nDiagEl
                exit
            endif
        enddo
    enddo
    allocate(diagElList(max(nDiagEl,1)), diagElF(3,8,6,max(nDiagEl,1)), diagElS(9,max(nDiagEl,1)))
    diagElList = 0
    diagElF = 0.0d0
    diagElS = 0.0d0
    do nel = 1, totalNumOfElements
        if (diagElIdx(nel) > 0) diagElList(diagElIdx(nel)) = nel
    enddo
    write(fname,'(A,I4.4,A)') 'diagEl.', me, '.txt'
    open(9934, file = trim(fname), status = 'replace')
    write(fname,'(A,I4.4,A)') 'diagElS.', me, '.txt'
    open(9935, file = trim(fname), status = 'replace')
    write(fname,'(A,I4.4,A)') 'diagNd.', me, '.txt'
    open(9936, file = trim(fname), status = 'replace')
    write(*,'(A,I4,A,I6,A,I6,A,3ES12.4,A,I3)') ' diagBudget rank ', me, ' pairs ', nDiagPr, ' elements ', nDiagEl, &
        ' rdampm rdampk dt ', rdampm, rdampk, dt, ' C_hg ', C_hg
end subroutine diagBudgetInit
'''

OUT_SUB = r'''
subroutine diagBudget(stage)
! stage 1: after assembly and the MPI exchange, before faulting (keep the
! assembled force). stage 2: after faulting, before the division by mass.
    use globalvar
    implicit none
    integer (kind = 4) :: stage, i, k, ie, ip, ift
    integer (kind = 8) :: nel, ns, nm, nd, es, em
    real (kind = dp) :: cen(3), fpS(3), fpM(3), fqS(3), fqM(3), msS, msM
    if (stage == 1) then
        diagFpre = nodalForceArr(1:totalNumOfEquations)
        return
    endif
    ift = 1
    do ie = 1, nDiagEl
        nel = diagElList(ie)
        cen = 0.0d0
        do k = 1, nen
            cen = cen + meshCoor(1:3,nodeElemIdRelation(k,nel))/8.0d0
        enddo
        write(9935,'(I7,ES24.15,I10,I4,12ES24.15)') nt, timeElapsed, nel, elemTypeArr(nel), cen, diagElS(1:9,ie)
        do k = 1, nen
            nd = nodeElemIdRelation(k,nel)
            write(9934,'(I7,ES24.15,I10,I4,3ES24.15,I3,3ES24.15,2I3,18ES24.15)') nt, timeElapsed, nel, elemTypeArr(nel), cen, k, &
                meshCoor(1:3,nd), diagNodeFlag(nd), numOfDofPerNodeArr(nd), diagElF(1:3,k,1:6,ie)
        enddo
    enddo
    do ip = 1, nDiagPr
        i = diagPrList(ip)
        ns = nsmp(1,i,ift)
        nm = nsmp(2,i,ift)
        do k = 1, 3
            es = eqNumIndexArr(eqNumStartIndexLoc(ns)+k)
            em = eqNumIndexArr(eqNumStartIndexLoc(nm)+k)
            fpS(k) = diagFpre(es)
            fpM(k) = diagFpre(em)
            fqS(k) = nodalForceArr(es)
            fqM(k) = nodalForceArr(em)
        enddo
        msS = nodalMassArr(es)
        msM = nodalMassArr(em)
        write(9936,'(I7,60ES24.15)') nt, timeElapsed, meshCoor(1:3,ns), arn(i,ift), fnms(ns), fnms(nm), msS, msM, &
            un(1:3,i,ift), us(1:3,i,ift), ud(1:3,i,ift), fpS, fpM, fqS, fqM, &
            velArr(1:3,ns), velArr(1:3,nm), dispArr(1:3,ns), dispArr(1:3,nm), &
            fric(FRIC_SLOT_TRACT_NORM:FRIC_SLOT_TRACT_NORM+2,i,ift), fric(FRIC_SLOT_INIT_NORM,i,ift), &
            fric(FRIC_SLOT_INIT_STRIKE_SHEAR,i,ift), fric(FRIC_SLOT_INIT_DIP_SHEAR,i,ift), dt
    enddo
    if (nt == nstep) then
        close(9934)
        close(9935)
        close(9936)
    endif
end subroutine diagBudget
'''
rep('driver.f90', "end subroutine driver\n", "end subroutine driver\n" + INIT_SUB + OUT_SUB)
print('edits applied')

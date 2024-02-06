! Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
! MIT
subroutine calcElemKU(shg,mate,vl,dl,stress,elresf,constk,porep,pstrmag,ex)
    use globalvar
    implicit none
    
    integer (kind = 4) :: i,j,j1,j2,j3
    real (kind = dp) :: constk,temp,porep
    real (kind = dp),dimension(nee) :: elresf,work,vl,dl
    real (kind = dp),dimension(nstr) :: strainrate,stressrate,strtemp,strain
    real (kind = dp)::stress(12),anestr(6),anestr1(6)
    real (kind = dp),dimension(nrowsh-1,nen) :: shg
    real (kind = dp),dimension(nrowb,nee) :: bb	!correspond to b
    real (kind = dp),dimension(nrowc,nrowc) :: cc	!correspond to c
    !...plasticity vlrables. B.D. 1/5/12
    real (kind = dp) :: strmea,taomax, yield, rjust,pstrmea,pstrmag,xc(3),ex(3,8)
    real (kind = dp),dimension(nstr) :: strdev,pstrinc
    !Parameters for Q model.
    real (kind = dp)::Qp,Qs,wkp,wks,taok,cv,cs
    real (kind = dp)::lam,miu,Mu,miuu,vols,mate(5)
    integer (kind = 4)::k,ip,iq,ir
    !
    pstrmag = 0.0d0
    stressrate = 0.0d0
    strain=0.0d0
    cc=0.0d0
    lam=mate(4)
    miu=mate(5)
    do i=1,3 
        cc(i,i)=lam+2*miu
        cc(i+3,i+3)=miu 
    enddo	
    cc(1,2)=lam 
    cc(2,1)=lam 
    cc(1,3)=lam
    cc(3,1)=lam 
    cc(2,3)=lam 
    cc(3,2)=lam 
    !...calcuate b from shg
    call qdcb(shg,bb)

    strainrate = 0.0d0	!initialize

    do i=1,nen
        j1 = ned * (i - 1) + 1
        j2 = ned * (i - 1) + 2
        j3 = ned * (i - 1) + 3      
        strainrate(1) = strainrate(1) + bb(1,j1) * vl(j1)
        strainrate(2) = strainrate(2) + bb(2,j2) * vl(j2)
        strainrate(3) = strainrate(3) + bb(3,j3) * vl(j3)
        strainrate(4) = strainrate(4) + bb(4,j2) * vl(j2) + bb(4,j3) * vl(j3)
        strainrate(5) = strainrate(5) + bb(5,j1) * vl(j1) + bb(5,j3) * vl(j3)
        strainrate(6) = strainrate(6) + bb(6,j1) * vl(j1) + bb(6,j2) * vl(j2)
        strain(1) = strain(1) + bb(1,j1) * dl(j1)
        strain(2) = strain(2) + bb(2,j2) * dl(j2)
        strain(3) = strain(3) + bb(3,j3) * dl(j3)
        strain(4) = strain(4) + bb(4,j2) * dl(j2) + bb(4,j3) * dl(j3)
        strain(5) = strain(5) + bb(5,j1) * dl(j1) + bb(5,j3) * dl(j3)
        strain(6) = strain(6) + bb(6,j1) * dl(j1) + bb(6,j2) * dl(j2)		
    enddo
    !...calculate stressrate
    ! Take into account zero in cc. B.D. 8/20/05
    do i=1,3
        do j=1,3
            stressrate(i) = stressrate(i) + cc(i,j)*strainrate(j)
        enddo
    enddo
    do i=4,6
        stressrate(i) = cc(i,i) * strainrate(i)
    enddo
    !...calculate stress from stress rate & previous stress. B.D. 1/5/12
    if (C_Q==0) then
        do i=1,6
            stress(i) = stress(i) + stressrate(i) * dt
            strdev(i) = stress(i)
        enddo
    elseif (C_Q==1) then
        xc=0.0d0	
        do i=1,3
            do j=1,8
                xc(i)=xc(i)+ex(i,j)
            enddo
        enddo
        xc=xc/8.0d0		
    !For DCPS		
        if (xc(3)>-1000.0d0)then
            Qs=10.0d0
            Qp=20.0d0
        else
            Qs=50.0d0
            Qp=100.0d0
        endif
    !For Tianjin
        ! if (mate(2)<1500)then
        ! Qs=0.02*mate(2)
        ! else
        ! Qs=0.1*mate(2)
        ! endif
        ! Qp=1.5*Qs
        ip=(xc(1)-(PMLb(2)+dx/2))/dx+1
        iq=(xc(2)-(PMLb(4)+dx/2))/dx+1
        ir=(xc(3)-(PMLb(5)+dx/2))/dx+1
        k=1+mod(ip,2)+2*mod(iq,2)+4*mod(ir,2)
        ! Modified Day and Bradley(2001) based on Liu(2006)	
        call qconstant(Qp,taok,wkp,k,cv)
        call qconstant(Qs,taok,wks,k,cs)
        wkp=wkp*8.0d0;
        wks=wks*8.0d0;
        miuu=miu*cs
        Mu=(lam+2*miu)*cv
        vols=strain(1)+strain(2)+strain(3)
        do i=1,6
            anestr1(i)=stress(i+6)
        enddo	
        do i=1,3
        anestr(i)=exp(-dt/taok)*anestr1(i)+(1-exp(-dt/taok)) &
             *(2*miuu*strain(i)*wks+(Mu*wkp-2*miuu*wks)*vols)
        enddo
        do i=4,6
        anestr(i)=exp(-dt/taok)*anestr1(i)+(1-exp(-dt/taok)) &
             *(miuu*strain(i)*wks)
        enddo	
        do i=1,6
            stress(i+6)=anestr(i)
        enddo
        stress(1)=2.0d0*miuu*strain(1)+(Mu-2.0d0*miuu)*vols-0.5d0*(anestr(1)+anestr1(1))
        stress(2)=2.0d0*miuu*strain(2)+(Mu-2.0d0*miuu)*vols-0.5d0*(anestr(2)+anestr1(2))	
        stress(3)=2.0d0*miuu*strain(3)+(Mu-2.0d0*miuu)*vols-0.5d0*(anestr(3)+anestr1(3))	
        stress(4)=2.0d0*miuu*strain(4)/2.0d0-0.5d0*(anestr(4)+anestr1(4))	
        stress(5)=2.0d0*miuu*strain(5)/2.0d0-0.5d0*(anestr(5)+anestr1(5))	
        stress(6)=2.0d0*miuu*strain(6)/2.0d0-0.5d0*(anestr(6)+anestr1(6))	
    endif!C_Q
    if (C_elastic==0) then
        !...Drucker-Prager plasticity in shear. B.D. 1/5/12
        ! refer to the benchmark problem description.
        strmea = (stress(1)+stress(2)+stress(3))/3.0d0
        do i=1,3
            strdev(i) = stress(i) - strmea
        enddo
        taomax = 0.5d0*(strdev(1)**2+strdev(2)**2+strdev(3)**2) &
            + strdev(4)**2+strdev(5)**2+strdev(6)**2  !second invlriant of deviator
        taomax=sqrt(taomax)  !kind of max shear
        yield=ccosphi-sinphi*(strmea + porep)  !yield stress
        if(yield<0.0d0) yield=0.0d0  !non-negative
        if(taomax > yield) then  !yielding, stress adjust in deviator domain
            !rjust = yield/taomax  !adjust ratio
            !implement viscoplasticity now. B.D. 6/2/12
            rjust=yield/taomax + (1-yield/taomax)*exp(-dt/tv)
            do i=1,6  !adjust stress
                stress(i) =strdev(i) * rjust
                !calculate plastic strain increment components
                pstrinc(i) = (strdev(i) - stress(i))/miu
                !back to stress domain by ading mean, which does not change
                if(i<=3) then
                    stress(i) = stress(i) + strmea
                endif
            enddo
            !calculate plastic strain increment scalar (amplitude)
            pstrmea = (pstrinc(1)+pstrinc(2)+pstrinc(3))/3.0d0
            do i=1,6
                pstrinc(i)=pstrinc(i) - pstrmea
            enddo
            pstrmag = 0.5d0*(pstrinc(1)**2+pstrinc(2)**2+pstrinc(3)**2) &
                    +pstrinc(4)**2+pstrinc(5)**2+pstrinc(6)**2
            pstrmag = sqrt(pstrmag)
        endif
     endif!C_elastic==0(Plastic)

    temp = constk * w
    do i=1,nstr
        strtemp(i) = temp * (stress(i) + rdampk * stressrate(i))
        !strtemp(i) = temp * stress(i)
    enddo

    do i=1,nen
        j1 = ned * (i - 1) + 1
        j2 = ned * (i - 1) + 2
        j3 = ned * (i - 1) + 3      
        work(j1) = bb(1,j1)*strtemp(1) + bb(5,j1)*strtemp(5) &
                + bb(6,j1)*strtemp(6)
        work(j2) = bb(2,j2)*strtemp(2) + bb(4,j2)*strtemp(4) &
                + bb(6,j2)*strtemp(6)
        work(j3) = bb(3,j3)*strtemp(3) + bb(4,j3)*strtemp(4) &
                + bb(5,j3)*strtemp(5)
    enddo		

    do i=1,nee
        elresf(i) = elresf(i) + work(i)
    enddo

end subroutine calcElemKU

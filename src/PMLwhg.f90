SUBROUTINE PMLwhg(vl,f,s,ex,mat1,shg,det,nel)
!efPML,evPML,esPML,ex,PMLb,c(1,1,1),eleshp(1,1,nel),det,dt
    use globalvar
    implicit none  
    !-------------------------------------------!
    ! Perfectly Matched Absorbing Boundary Layer.
    ! Input: 
    ! -vl: Velocity vector at step (n-1/2) on element nodes.
    ! -s: Stress vector at step (n) at the Guass Point of the element.   
    ! -f: Element force vector.
    ! -c: material matrix    
    ! -PMLb: PML boundary coordinates. 
    !     - PMLb(1):xmax2
    !     - PMLb(2):xmin2
    !     - PMLb(3):ymax2
    !    - PMLb(4):ymin2
    !    - PMLb(5):zmin2
    real(kind = dp),dimension(3,8)::vl    
    real(kind = dp),dimension(96)::f
    !vxx,vxy,vxz,
    !vyx,vyy,vyz,
    !vzx,vzy,vzz,
    !vhgx,vhgy,vhgz. *8 nodes.
    real(kind = dp),dimension(24)::va!total velocity
    real(kind = dp),dimension(21)::s
    !sxxx,sxxy,sxxz,
    !syyx,syyy,syyz,
    !szzx,szzy,szzz,
    !sxyx,sxyy,
    !sxzx,sxzz,
    !syzy,syzz.
    real(kind = dp)::lam,miu,det,vp,delta,kapa,rou,coef
    real(kind = dp)::Dx_vx,Dy_vy,Dz_vz,Dx_vy,Dy_vx,Dx_vz,Dz_vx,Dy_vz,Dz_vy
    real(kind = dp)::sxx,syy,szz,sxy,sxz,syz,s0(6)
    real(kind = dp),dimension(3,8)::ex
    real(kind = dp),dimension(3)::xc
    real(kind = dp),dimension(3)::damps
    real(kind = dp),dimension(3,8)::shg
    real(kind = dp),dimension(3,4)::q
    integer(kind=4),dimension(4,8)::fi
    real(kind = dp)::xmax2,xmin2,ymax2,ymin2,zmin2,&
                    maxdx,maxdy,maxdz,mat1(5)
    integer(kind=4)::i,j,k,j1,j2,j3,nel
    real (kind = dp),dimension(nrowb,nee) :: bb    !correspond to b
    real (kind = dp),dimension(nstr) :: strainrate,stressrate    
    real(kind = dp)::c(6,6)
    
    vp=mat1(1)
    lam=mat1(4)
    miu=mat1(5)
    !
    xmax2=PMLb(1)
    xmin2=PMLb(2)
    ymax2=PMLb(3)
    ymin2=PMLb(4)
    zmin2=PMLb(5)
    maxdx=PMLb(6)
    maxdy=PMLb(7)
    maxdz=PMLb(8)
    !
    xc=0.0d0
    do i=1,3
        do j=1,8
            xc(i)=xc(i)+ex(i,j)
        enddo
    enddo
    xc=xc/8

    ! Calculate damping profiles.
    if (xc(3)<zmin2) then !region 1
        damps(3)=abs(xc(3)-zmin2)        
        if (xc(1)>xmax2.and.xc(2)>ymax2) then !region 11
            damps(1)=abs(xc(1)-xmax2)
            damps(2)=abs(xc(2)-ymax2)            
        elseif (xc(1)>xmax2.and.xc(2)<ymin2) then !region 12
            damps(1)=abs(xc(1)-xmax2)
            damps(2)=abs(xc(2)-ymin2)    
        elseif (xc(1)<xmin2.and.xc(2)<ymin2) then !region 13
            damps(1)=abs(xc(1)-xmin2)
            damps(2)=abs(xc(2)-ymin2)            
        elseif (xc(1)<xmin2.and.xc(2)>xmax2) then !region 14
            damps(1)=abs(xc(1)-xmin2)
            damps(2)=abs(xc(2)-ymax2)
        elseif (xc(1)>xmax2.and.xc(2)>ymin2.and.xc(2)<ymax2) then !region 1_12
            damps(1)=abs(xc(1)-xmax2)
            damps(2)=0.0d0
        elseif (xc(2)<ymin2.and.xc(1)>xmin2.and.xc(1)<xmax2) then !region 1_23    
            damps(1)=0.0d0
            damps(2)=abs(xc(2)-ymin2)
        elseif (xc(1)<xmin2.and.xc(2)>ymin2.and.xc(2)<ymax2) then !region 1_34    
            damps(1)=abs(xc(1)-xmin2)
            damps(2)=0.0d0
        elseif (xc(2)>ymax2.and.xc(1)>xmin2.and.xc(1)<xmax2) then !region 1_41    
            damps(1)=0.0d0
            damps(2)=abs(xc(2)-ymax2)
        else
        !Middle area 9 missing previously.
        !Feb.18.2016/D.Liu
            damps(1)=0.0d0
            damps(2)=0.0d0
        endif
    elseif (xc(3)>zmin2) then !region 2
        damps(3)=0.0d0
        if (xc(1)>xmax2.and.xc(2)>ymax2) then !region 11
            damps(1)=abs(xc(1)-xmax2)
            damps(2)=abs(xc(2)-ymax2)            
        elseif (xc(1)>xmax2.and.xc(2)<ymin2) then !region 12
            damps(1)=abs(xc(1)-xmax2)
            damps(2)=abs(xc(2)-ymin2)    
        elseif (xc(1)<xmin2.and.xc(2)<ymin2) then !region 13
            damps(1)=abs(xc(1)-xmin2)
            damps(2)=abs(xc(2)-ymin2)            
        elseif (xc(1)<xmin2.and.xc(2)>xmax2) then !region 14
            damps(1)=abs(xc(1)-xmin2)
            damps(2)=abs(xc(2)-ymax2)
        elseif (xc(1)>xmax2.and.xc(2)>ymin2.and.xc(2)<ymax2) then !region 1_12
            damps(1)=abs(xc(1)-xmax2)
            damps(2)=0.0d0
        elseif (xc(2)<ymin2.and.xc(1)>xmin2.and.xc(1)<xmax2) then !region 1_23    
            damps(1)=0.0d0
            damps(2)=abs(xc(2)-ymin2)
        elseif (xc(1)<xmin2.and.xc(2)>ymin2.and.xc(2)<ymax2) then !region 1_34    
            damps(1)=abs(xc(1)-xmin2)
            damps(2)=0.0d0
        elseif (xc(2)>ymax2.and.xc(1)>xmin2.and.xc(1)<xmax2) then !region 1_41    
            damps(1)=0.0d0
            damps(2)=abs(xc(2)-ymax2)
        else
        !Middle area 9 missing previously.
        !Feb.18.2016/D.Liu
            damps(1)=0.0d0
            damps(2)=0.0d0
        endif
    endif
    do i=1,3
        if (i==1) then
        delta=nPML*maxdx
        elseif (i==2) then
        delta=nPML*maxdy
        elseif (i==3) then
        delta=nPML*maxdz        
        endif
        damps(i)=3*vmaxPML/2/delta*log(1/R)*(damps(i)/delta)**(2.)
    enddo    
    !-------------------------------------------------! 
    ! Calculate differentials of velocity.        
    do i=1,8
        va(3*(i-1)+1)=vl(1,i)
        va(3*(i-1)+2)=vl(2,i)
        va(3*(i-1)+3)=vl(3,i)
    enddo

    call qdcb(shg,bb)
    stressrate = 0.0d0
    strainrate = 0.0d0    !initialize
    do i=1,nen
        j1 = ned * (i - 1) + 1
        j2 = ned * (i - 1) + 2
        j3 = ned * (i - 1) + 3      
        strainrate(1) = strainrate(1) + bb(1,j1) * va(j1)
        strainrate(2) = strainrate(2) + bb(2,j2) * va(j2)
        strainrate(3) = strainrate(3) + bb(3,j3) * va(j3)
        strainrate(4) = strainrate(4) + bb(4,j2) * va(j2) + bb(4,j3) * va(j3)
        strainrate(5) = strainrate(5) + bb(5,j1) * va(j1) + bb(5,j3) * va(j3)
        strainrate(6) = strainrate(6) + bb(6,j1) * va(j1) + bb(6,j2) * va(j2)
    enddo
    c=0.0d0 
    do i=1,3 
        c(i,i)=lam+2.0d0*miu
        c(i+3,i+3)=miu 
    enddo    
    c(1,2)=lam 
    c(2,1)=lam 
    c(2,3)=lam 
    c(3,2)=lam 
    c(1,3)=lam 
    c(3,1)=lam 
    do i=1,3
        do j=1,3
            stressrate(i) = stressrate(i) + c(i,j)*strainrate(j)
        enddo
    enddo
    do i=4,6
        stressrate(i) = c(i,i) * strainrate(i)
    enddo    
    Dx_vx=0.0d0
    Dy_vy=0.0d0
    Dz_vz=0.0d0
    Dx_vy=0.0d0
    Dy_vx=0.0d0
    Dx_vz=0.0d0
    Dz_vx=0.0d0
    Dy_vz=0.0d0
    Dz_vy=0.0d0
    do i=1,8
        Dx_vx=Dx_vx+shg(1,i)*va(3*(i-1)+1)
        Dy_vy=Dy_vy+shg(2,i)*va(3*(i-1)+2)
        Dz_vz=Dz_vz+shg(3,i)*va(3*(i-1)+3)
        Dx_vy=Dx_vy+shg(1,i)*va(3*(i-1)+2)
        Dy_vx=Dy_vx+shg(2,i)*va(3*(i-1)+1)
        Dx_vz=Dx_vz+shg(1,i)*va(3*(i-1)+3)
        Dz_vx=Dz_vx+shg(3,i)*va(3*(i-1)+1)
        Dy_vz=Dy_vz+shg(2,i)*va(3*(i-1)+3)
        Dz_vy=Dz_vy+shg(3,i)*va(3*(i-1)+2)
    enddo
    !-------------------------------------------------!  
    ! Update stress
    s(1)=(lam+2*miu)*Dx_vx+(1/dt-damps(1)/2)*s(1)
    s(1)=s(1)/(1/dt+damps(1)/2)
    s(2)=lam*Dy_vy+(1/dt-damps(2)/2)*s(2)
    s(2)=s(2)/(1/dt+damps(2)/2)
    s(3)=lam*Dz_vz+(1/dt-damps(3)/2)*s(3)
    s(3)=s(3)/(1/dt+damps(3)/2)
    !-  
    s(4)=lam*Dx_vx+(1/dt-damps(1)/2)*s(4)
    s(4)=s(4)/(1/dt+damps(1)/2)
    s(5)=(lam+2*miu)*Dy_vy+(1/dt-damps(2)/2)*s(5)
    s(5)=s(5)/(1/dt+damps(2)/2)
    s(6)=lam*Dz_vz+(1/dt-damps(3)/2)*s(6)
    s(6)=s(6)/(1/dt+damps(3)/2)
    !-  
    s(7)=lam*Dx_vx+(1/dt-damps(1)/2)*s(7)
    s(7)=s(7)/(1/dt+damps(1)/2)
    s(8)=lam*Dy_vy+(1/dt-damps(2)/2)*s(8)
    s(8)=s(8)/(1/dt+damps(2)/2)
    s(9)=(lam+2*miu)*Dz_vz+(1/dt-damps(3)/2)*s(9)
    s(9)=s(9)/(1/dt+damps(3)/2)
    !-
    s(10)=miu*Dx_vy+(1/dt-damps(1)/2)*s(10)
    s(10)=s(10)/(1/dt+damps(1)/2)
     s(11)=miu*Dy_vx+(1/dt-damps(2)/2)*s(11)
    s(11)=s(11)/(1/dt+damps(2)/2) 
    !-    
    s(12)=miu*Dx_vz+(1/dt-damps(1)/2)*s(12)
    s(12)=s(12)/(1/dt+damps(1)/2)
     s(13)=miu*Dz_vx+(1/dt-damps(3)/2)*s(13)
    s(13)=s(13)/(1/dt+damps(3)/2) 
    !-    
    s(14)=miu*Dy_vz+(1/dt-damps(2)/2)*s(14)
    s(14)=s(14)/(1/dt+damps(2)/2)
     s(15)=miu*Dz_vy+(1/dt-damps(3)/2)*s(15)
    s(15)=s(15)/(1/dt+damps(3)/2)

    sxx=s(1)+s(2)+s(3)
    syy=s(4)+s(5)+s(6)    
    szz=s(7)+s(8)+s(9)
    sxy=s(10)+s(11)
    sxz=s(12)+s(13)
    syz=s(14)+s(15)    
    
    s0(1)=s(15+1)+rdampk * stressrate(1)
    s0(2)=s(15+2)+rdampk * stressrate(2)
    s0(3)=s(15+3)+rdampk * stressrate(3)
    s0(6)=s(15+6)+rdampk * stressrate(6)
    s0(5)=s(15+5)+rdampk * stressrate(5)
    s0(4)=s(15+4)+rdampk * stressrate(4)!!!

    !Calculate 'nodal forces'
    do i=1,8
        f((i-1)*12+1)=f((i-1)*12+1)-det*w*shg(1,i)*sxx
        f((i-1)*12+2)=f((i-1)*12+2)-det*w*shg(2,i)*sxy
        f((i-1)*12+3)=f((i-1)*12+3)-det*w*shg(3,i)*sxz
        f((i-1)*12+4)=f((i-1)*12+4)-det*w*shg(1,i)*sxy
        f((i-1)*12+5)=f((i-1)*12+5)-det*w*shg(2,i)*syy
        f((i-1)*12+6)=f((i-1)*12+6)-det*w*shg(3,i)*syz
        f((i-1)*12+7)=f((i-1)*12+7)-det*w*shg(1,i)*sxz
        f((i-1)*12+8)=f((i-1)*12+8)-det*w*shg(2,i)*syz
        f((i-1)*12+9)=f((i-1)*12+9)-det*w*shg(3,i)*szz
        f((i-1)*12+10)=f((i-1)*12+10)-&
        det*w*(bb(1,3*(i-1)+1)*s0(1)+bb(5,3*(i-1)+1)*s0(5)+bb(6,3*(i-1)+1)*s0(6))
        f((i-1)*12+11)=f((i-1)*12+11)-&
        det*w*(bb(2,3*(i-1)+2)*s0(2)+bb(4,3*(i-1)+2)*s0(4)+bb(6,3*(i-1)+2)*s0(6))
        f((i-1)*12+12)=f((i-1)*12+12)-&
        det*w*(bb(3,3*(i-1)+3)*s0(3)+bb(4,3*(i-1)+3)*s0(4)+bb(5,3*(i-1)+3)*s0(5))
    enddo
    !
end subroutine PMLwhg

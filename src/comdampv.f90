subroutine comdampv(x2,y,z,dv)
use globalvar
implicit none
integer(kind=4)::i
real (kind = dp) :: x2,y,z,xmax0,xmin0,ymax0,ymin0,zmin0,delta,maxdx,maxdy,maxdz
real (kind = dp),dimension(3)::damp
real (kind = dp),dimension(9)::dv
    !
    xmax0=PMLb(1)
    xmin0=PMLb(2)
    ymax0=PMLb(3)
    ymin0=PMLb(4)
    zmin0=PMLb(5)
    maxdx=PMLb(6)
    maxdy=PMLb(7)
    maxdz=PMLb(8)    
    !
    if (z<=zmin0) then !region 1
        damp(3)=abs(z-zmin0)        
        if (x2>=xmax0.and.y>=ymax0) then !region 11
            damp(1)=abs(x2-xmax0)
            damp(2)=abs(y-ymax0)            
        elseif (x2>=xmax0.and.y<=ymin0) then !region 12
            damp(1)=abs(x2-xmax0)
            damp(2)=abs(y-ymin0)    
        elseif (x2<=xmin0.and.y<=ymin0) then !region 13
            damp(1)=abs(x2-xmin0)
            damp(2)=abs(y-ymin0)            
        elseif (x2<=xmin0.and.y>=xmax0) then !region 14
            damp(1)=abs(x2-xmin0)
            damp(2)=abs(y-ymax0)
        elseif (x2>=xmax0.and.y>ymin0.and.y<ymax0) then !region 1_12
            damp(1)=abs(x2-xmax0)
            damp(2)=0.0d0
        elseif (y<=ymin0.and.x2>xmin0.and.x2<xmax0) then !region 1_23    
            damp(1)=0.0d0
            damp(2)=abs(y-ymin0)
        elseif (x2<=xmin0.and.y>ymin0.and.y<ymax0) then !region 1_34    
            damp(1)=abs(x2-xmin0)
            damp(2)=0.0d0
        elseif (y>=ymax0.and.x2>xmin0.and.x2<xmax0) then !region 1_41    
            damp(1)=0.0d0
            damp(2)=abs(y-ymax0)
        else
        !Middle area 9 missing previously.
        !Feb.18.2016/D.Liu
            damp(1)=0.0d0
            damp(2)=0.0d0 
        endif
    elseif (z>zmin0) then !region 2
        damp(3)=0.0d0
        if (x2>=xmax0.and.y>=ymax0) then !region 11
            damp(1)=abs(x2-xmax0)
            damp(2)=abs(y-ymax0)            
        elseif (x2>=xmax0.and.y<=ymin0) then !region 12
            damp(1)=abs(x2-xmax0)
            damp(2)=abs(y-ymin0)    
        elseif (x2<=xmin0.and.y<=ymin0) then !region 13
            damp(1)=abs(x2-xmin0)
            damp(2)=abs(y-ymin0)            
        elseif (x2<=xmin0.and.y>=xmax0) then !region 14
            damp(1)=abs(x2-xmin0)
            damp(2)=abs(y-ymax0)
        elseif (x2>=xmax0.and.y>ymin0.and.y<ymax0) then !region 1_12
            damp(1)=abs(x2-xmax0)
            damp(2)=0.0d0
        elseif (y<=ymin0.and.x2>xmin0.and.x2<xmax0) then !region 1_23    
            damp(1)=0.0d0
            damp(2)=abs(y-ymin0)
        elseif (x2<=xmin0.and.y>ymin0.and.y<ymax0) then !region 1_34    
            damp(1)=abs(x2-xmin0)
            damp(2)=0.0d0
        elseif (y>=ymax0.and.x2>xmin0.and.x2<xmax0) then !region 1_41    
            damp(1)=0.0d0
            damp(2)=abs(y-ymax0)
        else 
        !Middle area 9 missing previously.
        !Feb.18.2016/D.Liu
        !Actually this region does not exist
            damp(1)=0.0d0 
            damp(2)=0.0d0 
        endif
    endif    
!For Double-couple point source. Ma and Liu (2006)     
    ! if (z>-1000.) then
    ! vp=2800.
    ! elseif(z==-1000.)then
    ! vp=4400.
    ! else
    ! vp=6000.
    ! endif
!For TPV 8
!    vp=5716.
!For TianJin
    do i=1,3
        if (i==1) then
        delta=nPML*maxdx
        elseif (i==2) then
        delta=nPML*maxdy
        elseif (i==3) then
        delta=nPML*maxdz        
        endif    
        damp(i)=3.0d0*vmaxPML/2.0d0/delta*log(1.0d0/R)*(damp(i)/delta)**(2.0d0)
    enddo
    dv(1)=damp(1)
    dv(2)=damp(2)
    dv(3)=damp(3)
    dv(4)=damp(1)
    dv(5)=damp(2)
    dv(6)=damp(3)
    dv(7)=damp(1)
    dv(8)=damp(2)
    dv(9)=damp(3)
    do i=1,9
        if (dv(i)<0.0d0) then
            write(*,*) 'wrong dv'
            stop
        endif
    enddo
end subroutine comdampv

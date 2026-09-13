subroutine computePMLDampingVector(x2,y,z,dv)
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
    call pmlRegionDistance(x2, y, z, xmax0, xmin0, ymax0, ymin0, zmin0, .true., damp)
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
end subroutine computePMLDampingVector

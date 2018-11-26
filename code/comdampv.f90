subroutine comdampv(x,y,z,PMLb,dv)
use globalvar
implicit none
integer(kind=4)::i
real (kind=8) :: x,y,z,xmax,xmin,ymax,ymin,zmin,vp,delta,maxdx,maxdy,maxdz
real (kind=8),dimension(8):: PMLb
real (kind=8),dimension(3)::damp
real (kind=8),dimension(9)::dv
	!
	xmax=PMLb(1)
	xmin=PMLb(2)
	ymax=PMLb(3)
	ymin=PMLb(4)
	zmin=PMLb(5)
	maxdx=PMLb(6)
	maxdy=PMLb(7)
	maxdz=PMLb(8)
	!
	if (z<=zmin) then !region 1
		damp(3)=abs(z-zmin)		
		if (x>=xmax.and.y>=ymax) then !region 11
			damp(1)=abs(x-xmax)
			damp(2)=abs(y-ymax)			
		elseif (x>=xmax.and.y<=ymin) then !region 12
			damp(1)=abs(x-xmax)
			damp(2)=abs(y-ymin)	
		elseif (x<=xmin.and.y<=ymin) then !region 13
			damp(1)=abs(x-xmin)
			damp(2)=abs(y-ymin)			
		elseif (x<=xmin.and.y>=xmax) then !region 14
			damp(1)=abs(x-xmin)
			damp(2)=abs(y-ymax)
		elseif (x>=xmax.and.y>ymin.and.y<ymax) then !region 1_12
			damp(1)=abs(x-xmax)
			damp(2)=0.0
		elseif (y<=ymin.and.x>xmin.and.x<xmax) then !region 1_23	
			damp(1)=0.0
			damp(2)=abs(y-ymin)
		elseif (x<=xmin.and.y>ymin.and.y<ymax) then !region 1_34	
			damp(1)=abs(x-xmin)
			damp(2)=0.0
		elseif (y>=ymax.and.x>xmin.and.x<xmax) then !region 1_41	
			damp(1)=0.0
			damp(2)=abs(y-ymax)
		endif
	elseif (z>zmin) then !region 2
		damp(3)=0.0
		if (x>=xmax.and.y>=ymax) then !region 11
			damp(1)=abs(x-xmax)
			damp(2)=abs(y-ymax)			
		elseif (x>=xmax.and.y<=ymin) then !region 12
			damp(1)=abs(x-xmax)
			damp(2)=abs(y-ymin)	
		elseif (x<=xmin.and.y<=ymin) then !region 13
			damp(1)=abs(x-xmin)
			damp(2)=abs(y-ymin)			
		elseif (x<=xmin.and.y>=xmax) then !region 14
			damp(1)=abs(x-xmin)
			damp(2)=abs(y-ymax)
		elseif (x>=xmax.and.y>ymin.and.y<ymax) then !region 1_12
			damp(1)=abs(x-xmax)
			damp(2)=0.0
		elseif (y<=ymin.and.x>xmin.and.x<xmax) then !region 1_23	
			damp(1)=0.0
			damp(2)=abs(y-ymin)
		elseif (x<=xmin.and.y>ymin.and.y<ymax) then !region 1_34	
			damp(1)=abs(x-xmin)
			damp(2)=0.0
		elseif (y>=ymax.and.x>xmin.and.x<xmax) then !region 1_41	
			damp(1)=0.0
			damp(2)=abs(y-ymax)
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
	vp=5716.
	do i=1,3
		if (i==1) then
		delta=nPML*maxdx
		elseif (i==2) then
		delta=nPML*maxdy
		elseif (i==3) then
		delta=nPML*maxdz		
		endif	
		damp(i)=3*vp/2/delta*log(1/R)*(damp(i)/delta)**(2.)
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
		if (dv(i)<0.0) then
			write(*,*) 'wrong dv'
			stop
		endif
	enddo
end subroutine comdampv

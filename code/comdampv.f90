subroutine comdampv(x2,y,z,dv)
use globalvar
implicit none
integer(kind=4)::i
real (kind=8) :: x2,y,z,xmax1,xmin1,ymax1,ymin1,zmin1,delta,maxdx,maxdy,maxdz
real (kind=8),dimension(3)::damp
real (kind=8),dimension(9)::dv
	!
	xmax1=PMLb(1)
	xmin1=PMLb(2)
	ymax1=PMLb(3)
	ymin1=PMLb(4)
	zmin1=PMLb(5)
	maxdx=PMLb(6)
	maxdy=PMLb(7)
	maxdz=PMLb(8)	
	!
	if (z<=zmin1) then !region 1
		damp(3)=abs(z-zmin1)		
		if (x2>=xmax1.and.y>=ymax1) then !region 11
			damp(1)=abs(x2-xmax1)
			damp(2)=abs(y-ymax1)			
		elseif (x2>=xmax1.and.y<=ymin1) then !region 12
			damp(1)=abs(x2-xmax1)
			damp(2)=abs(y-ymin1)	
		elseif (x2<=xmin1.and.y<=ymin1) then !region 13
			damp(1)=abs(x2-xmin1)
			damp(2)=abs(y-ymin1)			
		elseif (x2<=xmin1.and.y>=xmax1) then !region 14
			damp(1)=abs(x2-xmin1)
			damp(2)=abs(y-ymax1)
		elseif (x2>=xmax1.and.y>ymin1.and.y<ymax1) then !region 1_12
			damp(1)=abs(x2-xmax1)
			damp(2)=0.0d0
		elseif (y<=ymin1.and.x2>xmin1.and.x2<xmax1) then !region 1_23	
			damp(1)=0.0d0
			damp(2)=abs(y-ymin1)
		elseif (x2<=xmin1.and.y>ymin1.and.y<ymax1) then !region 1_34	
			damp(1)=abs(x2-xmin1)
			damp(2)=0.0d0
		elseif (y>=ymax1.and.x2>xmin1.and.x2<xmax1) then !region 1_41	
			damp(1)=0.0d0
			damp(2)=abs(y-ymax1)
		else
		!Middle area 9 missing previously.
		!Feb.18.2016/D.Liu
			damp(1)=0.0d0
			damp(2)=0.0d0 
		endif
	elseif (z>zmin1) then !region 2
		damp(3)=0.0d0
		if (x2>=xmax1.and.y>=ymax1) then !region 11
			damp(1)=abs(x2-xmax1)
			damp(2)=abs(y-ymax1)			
		elseif (x2>=xmax1.and.y<=ymin1) then !region 12
			damp(1)=abs(x2-xmax1)
			damp(2)=abs(y-ymin1)	
		elseif (x2<=xmin1.and.y<=ymin1) then !region 13
			damp(1)=abs(x2-xmin1)
			damp(2)=abs(y-ymin1)			
		elseif (x2<=xmin1.and.y>=xmax1) then !region 14
			damp(1)=abs(x2-xmin1)
			damp(2)=abs(y-ymax1)
		elseif (x2>=xmax1.and.y>ymin1.and.y<ymax1) then !region 1_12
			damp(1)=abs(x2-xmax1)
			damp(2)=0.0d0
		elseif (y<=ymin1.and.x2>xmin1.and.x2<xmax1) then !region 1_23	
			damp(1)=0.0d0
			damp(2)=abs(y-ymin1)
		elseif (x2<=xmin1.and.y>ymin1.and.y<ymax1) then !region 1_34	
			damp(1)=abs(x2-xmin1)
			damp(2)=0.0d0
		elseif (y>=ymax1.and.x2>xmin1.and.x2<xmax1) then !region 1_41	
			damp(1)=0.0d0
			damp(2)=abs(y-ymax1)
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
!	vp=5716.
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

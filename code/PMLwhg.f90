SUBROUTINE PMLwhg(va,f,v,s,ex,PMLb,c,shg,det,nt,bdamp,nel,me)!efPML,evPML,esPML,ex,PMLb,c(1,1,1),eleshp(1,1,nel),det,dt
	use globalvar
	implicit none  
	! D.L. Feb/2015
	!-------------------------------------------!
	! Perfectly Matched Absorbing Boundary Layer.
	! Input: 
	! -vl: Total velocity vector at step (n-1/2) on element nodes.
	! -s: Stress vector at step (n) at the Guass Point of the element.   
	! -f: Element force vector.
	! -dt: Time interval
	! -c: Elastic material matrix	
	! -PMLb: PML boundary coordinates. 
	! 	- PMLb(1):xmax
	! 	- PMLb(2):xmin
	! 	- PMLb(3):ymax
	!	- PMLb(4):ymin
	!	- PMLb(5):zmin
	! 	- PMLb(6):hx
	!	- PMLb(7):hy
	!	- PMLb(8):hz
	real(kind=8),dimension(96)::f,v
	real(kind=8),dimension(24)::va!Total velocity
	!vxx,vxy,vxz,
	!vyx,vyy,vyz,
	!vzx,vzy,vzz,
	!vhgx,vhgy,vhgz. *8 nodes.
	real(kind=8),dimension(21)::s
	!sxxx,sxxy,sxxz,
	!syyx,syyy,syyz,
	!szzx,szzy,szzz,
	!sxyx,sxyy,
	!sxzx,sxzz,
	!syzy,syzz.
	real(kind=8)::lam,miu,det,R,vp,delta,kapa,rou,coef
	real(kind=8)::Dx_vx,Dy_vy,Dz_vz,Dx_vy,Dy_vx,Dx_vz,Dz_vx,Dy_vz,Dz_vy
	real(kind=8)::sxx,syy,szz,sxy,sxz,syz,s0(6)
	real(kind=8),dimension(6,6)::c
	real(kind=8),dimension(3,8)::ex
	real(kind=8),dimension(3)::xc
	real(kind=8),dimension(3)::damps
	real(kind=8),dimension(3,8)::shg
	real(kind=8),dimension(3,4)::q
	integer(kind=4),dimension(4,8)::fi
	real(kind=8)::xmax,xmin,ymax,ymin,zmin,PMLb(8),maxdx,maxdy,maxdz,bdamp
	integer(kind=4)::i,j,k,nt,j1,j2,j3,nel,me
	real (kind=8),dimension(nrowb,nee) :: bb
	real (kind=8),dimension(nstr) :: strainrate,stressrate	
	lam=c(1,2)
	miu=c(4,4)
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
	xc=0.0
	do i=1,3
		do j=1,8
			xc(i)=xc(i)+ex(i,j)
		enddo
	enddo
	xc=xc/8
	! Calculate damping profiles.
	if (xc(3)<zmin) then !region 1
		damps(3)=abs(xc(3)-zmin)		
		if (xc(1)>xmax.and.xc(2)>ymax) then !region 11
			damps(1)=abs(xc(1)-xmax)
			damps(2)=abs(xc(2)-ymax)			
		elseif (xc(1)>xmax.and.xc(2)<ymin) then !region 12
			damps(1)=abs(xc(1)-xmax)
			damps(2)=abs(xc(2)-ymin)	
		elseif (xc(1)<xmin.and.xc(2)<ymin) then !region 13
			damps(1)=abs(xc(1)-xmin)
			damps(2)=abs(xc(2)-ymin)			
		elseif (xc(1)<xmin.and.xc(2)>xmax) then !region 14
			damps(1)=abs(xc(1)-xmin)
			damps(2)=abs(xc(2)-ymax)
		elseif (xc(1)>xmax.and.xc(2)>ymin.and.xc(2)<ymax) then !region 1_12
			damps(1)=abs(xc(1)-xmax)
			damps(2)=0.0
		elseif (xc(2)<ymin.and.xc(1)>xmin.and.xc(1)<xmax) then !region 1_23	
			damps(1)=0.0
			damps(2)=abs(xc(2)-ymin)
		elseif (xc(1)<xmin.and.xc(2)>ymin.and.xc(2)<ymax) then !region 1_34	
			damps(1)=abs(xc(1)-xmin)
			damps(2)=0.0
		elseif (xc(2)>ymax.and.xc(1)>xmin.and.xc(1)<xmax) then !region 1_41	
			damps(1)=0.0
			damps(2)=abs(xc(2)-ymax)
		endif
	elseif (xc(3)>zmin) then !region 2
		damps(3)=0.0
		if (xc(1)>xmax.and.xc(2)>ymax) then !region 11
			damps(1)=abs(xc(1)-xmax)
			damps(2)=abs(xc(2)-ymax)			
		elseif (xc(1)>xmax.and.xc(2)<ymin) then !region 12
			damps(1)=abs(xc(1)-xmax)
			damps(2)=abs(xc(2)-ymin)	
		elseif (xc(1)<xmin.and.xc(2)<ymin) then !region 13
			damps(1)=abs(xc(1)-xmin)
			damps(2)=abs(xc(2)-ymin)			
		elseif (xc(1)<xmin.and.xc(2)>xmax) then !region 14
			damps(1)=abs(xc(1)-xmin)
			damps(2)=abs(xc(2)-ymax)
		elseif (xc(1)>xmax.and.xc(2)>ymin.and.xc(2)<ymax) then !region 1_12
			damps(1)=abs(xc(1)-xmax)
			damps(2)=0.0
		elseif (xc(2)<ymin.and.xc(1)>xmin.and.xc(1)<xmax) then !region 1_23	
			damps(1)=0.0
			damps(2)=abs(xc(2)-ymin)
		elseif (xc(1)<xmin.and.xc(2)>ymin.and.xc(2)<ymax) then !region 1_34	
			damps(1)=abs(xc(1)-xmin)
			damps(2)=0.0
		elseif (xc(2)>ymax.and.xc(1)>xmin.and.xc(1)<xmax) then !region 1_41	
			damps(1)=0.0
			damps(2)=abs(xc(2)-ymax)
		endif
	endif
	R=0.01
	vp=6000.
	do i=1,3
		if (i==1) then
			delta=6*maxdx
		elseif (i==2) then
			delta=6*maxdy
		elseif (i==3) then
			delta=6*maxdz		
		endif
		damps(i)=3*vp/2/delta*log(1/R)*(damps(i)/delta)**(2.)
	enddo	
	!-------------------------------------------------! 
	! Calculate differentials of velocity.		
	stressrate = 0.0
	!...calcuate b from shg
	call qdcb(shg,bb)
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
	do i=1,3
		do j=1,3
			stressrate(i) = stressrate(i) + c(i,j)*strainrate(j)
		enddo
	enddo
	do i=4,6
		stressrate(i) = c(i,i) * strainrate(i)
	enddo	
	Dx_vx=0.0
	Dy_vy=0.0
	Dz_vz=0.0
	Dx_vy=0.0
	Dy_vx=0.0
	Dx_vz=0.0
	Dz_vx=0.0
	Dy_vz=0.0
	Dz_vy=0.0
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
	! Take account of damping propertional to K.  
	s0(1)=s(15+1)+bdamp * stressrate(1)
	s0(2)=s(15+2)+bdamp * stressrate(2)
	s0(3)=s(15+3)+bdamp * stressrate(3)
	s0(4)=s(15+4)+bdamp * stressrate(4)
	s0(5)=s(15+5)+bdamp * stressrate(5)	
	s0(6)=s(15+6)+bdamp * stressrate(6)
	!Calculate nodal forces.
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
		! Nodal force contribution from initial stress. 
		f((i-1)*12+10)=f((i-1)*12+10)-det*w*(bb(1,3*(i-1)+1)*s0(1)+bb(5,3*(i-1)+1)*s0(5)+bb(6,3*(i-1)+1)*s0(6))
		f((i-1)*12+11)=f((i-1)*12+11)-det*w*(bb(2,3*(i-1)+2)*s0(2)+bb(4,3*(i-1)+2)*s0(4)+bb(6,3*(i-1)+2)*s0(6))
		f((i-1)*12+12)=f((i-1)*12+12)-det*w*(bb(3,3*(i-1)+3)*s0(3)+bb(4,3*(i-1)+3)*s0(4)+bb(5,3*(i-1)+3)*s0(5))
	enddo
	!if(xc(1)<51.and.xc(1)>49.and.xc(2)<51.and.xc(2)>49.and.(xc(3)>(zmin-200))) then
	!write(*,*) 'x,y,z',xc(1),xc(2),xc(3)
	!write(*,*) 'nel==',nel,me
	!write(*,*) 'PMLvel1',va(1),va(2),va(3)
	!write(*,*) 'PMLvel2',va(4),va(5),va(6)
	!write(*,*) 'PMLvel3',va(7),va(8),va(9)
	!write(*,*) 'PMLvel4',va(10),va(11),va(12)
	!write(*,*) 'PMLvel5',va(13),va(14),va(15)
	!write(*,*) 'PMLvel6',va(16),va(17),va(18)
	!write(*,*) 'PMLvel7',va(19),va(20),va(21)
	!write(*,*) 'PMLvel8',va(22),va(23),va(24)
	!write(*,*) 'PMLstr',sxx,syy,szz,sxy,sxz,syz
	!endif	
	!-------------------------------------------------! 
	!viscous hourglass control
	kapa=0.05!0.05~0.15
	rou=2670.
	coef=kapa*rou*vp*(det*w)**(2./3.)
	fi(1,1)=1;fi(1,2)=1;fi(1,3)=-1;fi(1,4)=-1;fi(1,5)=-1;fi(1,6)=-1;fi(1,7)=1;fi(1,8)=1
	fi(2,1)=1;fi(2,2)=-1;fi(2,3)=-1;fi(2,4)=1;fi(2,5)=-1;fi(2,6)=1;fi(2,7)=1;fi(2,8)=-1
	fi(3,1)=1;fi(3,2)=-1;fi(3,3)=1;fi(3,4)=-1;fi(3,5)=1;fi(3,6)=-1;fi(3,7)=1;fi(3,8)=-1
	fi(4,1)=1;fi(4,2)=-1;fi(4,3)=1;fi(4,4)=-1;fi(4,5)=-1;fi(4,6)=1;fi(4,7)=-1;fi(4,8)=1
	!fi1=[1,1,-1,-1,-1,-1,1,1]
	!fi2=[1,-1,-1,1,-1,1,1,-1]
	!fi3=[1,-1,1,-1,1,-1,1,-1]
	!fi4=[1,-1,1,-1,-1,1,-1,1]
	q=0.0!q(3,4)
	do i=1,3
		do j=1,4
			do k=1,8
				!q(i,j)=q(i,j)+va(3*(k-1)+i)*fi(j,k)
			enddo
		enddo
	enddo
	do k=1,8
		do i=1,3
			do j=1,4
				!f((k-1)*12+9+i)=f((k-1)*12+9+i)-coef*q(i,j)*fi(j,k)
			enddo
		enddo
	enddo	
end subroutine PMLwhg
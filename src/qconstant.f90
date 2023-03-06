subroutine qconstant(Q,rtaok,rwk,k,c1)
implicit none
real(kind = dp)::taok(8),alfk(8),betk(8)
real(kind = dp)::kapa,Q,rtaok,rwk,a1,b1,c1,ref,ak0,bk0
integer(kind=4)::k,i
taok=(/1.72333e-3,1.80701e-3,5.38887e-3,1.99322e-2,8.49833e-2,4.09335e-1,2.05951,13.2629/)
alfk=(/1.66958e-2,3.81644e-2,9.84666e-3,-1.36803e-2,-2.85125e-2,-5.37309e-2,-6.65035e-2,-1.33696e-1/) 
betk=(/8.98758e-2,6.84635e-2,9.67052e-2,1.20172e-1,1.30728e-1,1.38746e-1,1.40705e-1,2.14647e-1/) 
kapa=3.071+1.433*Q**(-1.158)*log(Q/5)
kapa=kapa/(1+0.415*Q)!EQ 6b for Q from 5 to 5000 (Liu and Archuleta,2006)
!Comments: frequency band 0.01~50Hz
rwk=kapa*(kapa*alfk(k)+betk(k))
rtaok=taok(k)
a1=1.0
b1=0.0
ref=2*3.1415926
! do i=1,8
	! a1=a1-kapa*(kapa*alfk(i)+betk(i))/(1+(ref*taok(i))**2)
	! b1=b1+kapa*(kapa*alfk(i)+betk(i))*ref*taok(i)/(1+(ref*taok(i))**2)
! enddo
! c1=sqrt(a1**2+b1**2)
ak0=1-rwk*8/(1+(rtaok*ref)**2)
bk0=rwk*8*ref*rtaok/(1+(rtaok*ref)**2)
c1=0.5*(ak0**2+bk0**2)**(-0.5)
c1=c1*(1+ak0*(ak0**2+bk0**2)**(-0.5))
! if (c1>=1.0) then
! write(*,*) 'wrong c1'
	! stop 555
! endif
end subroutine qconstant
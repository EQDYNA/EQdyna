SUBROUTINE contm(shg,det,elmass,constm)
use globalvar
implicit none     
!=============================================================
!### program to form mass matrix for a continuum element
!        with "nen" nodes: lumped mass only here!
!=============================================================
integer (kind=4) :: j,n,k
real (kind=8) :: det,dsum,totmas,constm,temp1,temp2
real (kind=8),dimension(nee) :: elmass
real (kind=8),dimension(nen) :: work
real (kind=8),dimension(nrowsh,nen) :: shg
!...initialize
dsum   = 0.0
totmas = 0.0
work   = 0.0
!...calculate
totmas = constm*w*det
do j=1,nen
	temp2 = totmas*shg(nrowsh,j)**2
    dsum = dsum + temp2
    work(j) = work(j) + temp2
enddo
!...scale diagonal to conserve total mass
temp1 = totmas/dsum
!...store terms in a column
do j=1,nen
	temp2 = temp1*work(j)
    n = (j - 1)*ned
    do k=1,ned
		elmass(n + k) = temp2
		!it seems 333 components are same for one node.
		!B.D. 2/23/05
    enddo
enddo
!
end SUBROUTINE contm
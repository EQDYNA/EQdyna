SUBROUTINE qdcshg(xl,det,shl,shg,nel,xs,lcubic)
  use globalvar
  implicit none
  !
  !...program to calculate global derivatives of shape function and
  !	Jacobian determinat for a 8-node hexahedral element.
  !	1 Gaussian point only!
  !
  logical :: lcubic
  integer(kind=4) :: i,j,k,nel
  real(kind=8) :: det,temp,cof11,cof12,cof13,cof21,cof22,cof23, &
		cof31,cof32,cof33
  real(kind=8),dimension(3,3) :: xs
  real(kind=8),dimension(nrowsh,nen) :: shl,shg,shlg
  real(kind=8),dimension(nesd,nen) :: xl
  !
  !...equal local to global ,first
  shg = shl
  !...deal with wedge degeneration. B.D. 11/27/08
  if(.not.lcubic) then	!degeneration to wedge element by following the book
    do i=1,nrowsh
      shg(i,3) = shl(i,3) + shl(i,4)	!always 4=3, 8=7
      shg(i,4) = 0.0
      shg(i,7) = shl(i,7) + shl(i,8)
      shg(i,8) = 0.0
    enddo
  endif
  !...calculate x,s
  do i=1,3
    do j=1,3
      temp = 0.0
      do k=1,nen
      temp = temp + shg(i,k) * xl(j,k)
      enddo
      xs(j,i) = temp
    enddo
  enddo
  !...cofactors (see book p149-150)
  cof11 = xs(2,2)*xs(3,3) - xs(2,3)*xs(3,2)
  cof12 = xs(2,3)*xs(3,1) - xs(2,1)*xs(3,3)
  cof13 = xs(2,1)*xs(3,2) - xs(2,2)*xs(3,1)
  cof21 = xs(3,2)*xs(1,3) - xs(3,3)*xs(1,2)
  cof22 = xs(3,3)*xs(1,1) - xs(3,1)*xs(1,3)
  cof23 = xs(3,1)*xs(1,2) - xs(3,2)*xs(1,1)
  cof31 = xs(1,2)*xs(2,3) - xs(1,3)*xs(2,2)
  cof32 = xs(1,3)*xs(2,1) - xs(1,1)*xs(2,3)
  cof33 = xs(1,1)*xs(2,2) - xs(1,2)*xs(2,1)
  !...determinant
  det = xs(1,1)*cof11 + xs(1,2)*cof12 + xs(1,3)*cof13
  if(det <= 0.0) then
    write(*,1000) nel
    write(*,*) 'det=',det
    stop	!non-positive det: terminate
  endif
  !...derivatives of global shape function
  shlg = shg
  do i=1,nen
    shg(1,i) = (shlg(1,i)*cof11 + shlg(2,i)*cof12  &
		+ shlg(3,i)*cof13)/det
    shg(2,i) = (shlg(1,i)*cof21 + shlg(2,i)*cof22  &
		+ shlg(3,i)*cof23)/det
    shg(3,i) = (shlg(1,i)*cof31 + shlg(2,i)*cof32  &
		+ shlg(3,i)*cof33)/det
  enddo
  !...now,get s,x from cof and save in xs to return
  xs=reshape((/cof11,cof12,cof13, cof21,cof22,cof23, &
  	cof31,cof32,cof33/),(/3,3/))
  xs = xs / det	
!
1000 format('1','non-positive determinant in element number  ',i10)
!
end SUBROUTINE qdcshg

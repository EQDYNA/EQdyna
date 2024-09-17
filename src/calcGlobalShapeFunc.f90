subroutine calcGlobalShapeFunc(xl,det,globalShapeFunc,nel,xs,lcubic)
  use globalvar
  implicit none
  ! calculate global derivatives of shape function and
  ! Jacobian determinant for a reduced order hexahedral element.
  logical :: lcubic
  
  integer (kind = 4) :: i,j,k,nel
  
  real(kind = dp) :: det,temp, &
                    cof11,cof12,cof13, &
                    cof21,cof22,cof23, &
                    cof31,cof32,cof33
  real(kind = dp), dimension(3,3)        :: xs
  real(kind = dp), dimension(nrowsh,nen) :: globalShapeFunc,tmpGlobalShapeFunc
  real(kind = dp), dimension(nesd,nen)   :: xl

  !...equal local to global ,first
  globalShapeFunc = localShapeFunc
  ! degeneration; Hughes, P125, 
  if( elemTypeArr(nel)==11 .or. elemTypeArr(nel)==12 ) then    !degeneration to wedge element by following the book
    do i=1,nrowsh
      globalShapeFunc(i,3) = localShapeFunc(i,3) + localShapeFunc(i,4)  
      globalShapeFunc(i,4) = 0.0d0
      globalShapeFunc(i,7) = localShapeFunc(i,7) + localShapeFunc(i,8)
      globalShapeFunc(i,8) = 0.0d0
    enddo
  endif
  !...calculate x,s
  do i=1,3
    do j=1,3
      temp = 0.0d0
      do k=1,nen
      temp = temp + globalShapeFunc(i,k) * xl(j,k)
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
  ! Calculate the determinant
  det   = xs(1,1)*cof11 + xs(1,2)*cof12 + xs(1,3)*cof13
  
  if (det <= 0.0d0) then
    write(*,*) 'Non-positive determinant; stop the code.'
    write(*,*) 'Element id is ', nel
    write(*,*) 'det=', det
    stop
  endif
  
  !Calculate derivatives of global shape function
  tmpGlobalShapeFunc  = globalShapeFunc
  
  do i  = 1,nen
    globalShapeFunc(1,i) = (tmpGlobalShapeFunc(1,i)*cof11 + tmpGlobalShapeFunc(2,i)*cof12  &
                + tmpGlobalShapeFunc(3,i)*cof13)/det
    globalShapeFunc(2,i) = (tmpGlobalShapeFunc(1,i)*cof21 + tmpGlobalShapeFunc(2,i)*cof22  &
                + tmpGlobalShapeFunc(3,i)*cof23)/det
    globalShapeFunc(3,i) = (tmpGlobalShapeFunc(1,i)*cof31 + tmpGlobalShapeFunc(2,i)*cof32  &
                + tmpGlobalShapeFunc(3,i)*cof33)/det
  enddo
  !...now,get s,x from cof and save in xs to return
  xs    = reshape((/cof11,cof12,cof13, cof21,cof22,cof23, &
                 cof31,cof32,cof33/),(/3,3/))
  xs    = xs/det

end subroutine calcGlobalShapeFunc

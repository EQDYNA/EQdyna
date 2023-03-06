SUBROUTINE contma(elmass,al,elresf)
  use globalvar
  implicit none
  !
  !### program to calculate inertial and gravity/body force ("-m*(a-g)")
  !	for a continuum element with "nen" nodes. In this version,
  !	for lumped mass only.
  !
  integer (kind=4) :: i,j,k
  real (kind = dp),dimension(nee) :: elmass
  real (kind = dp),dimension(ned,nen) :: elresf,al
  do j=1,nen
    k = (j - 1)*ned
    do i=1,ned
      elresf(i,j) = elresf(i,j) - al(i,j)*elmass(k + i)
      !note: - m * a is accumulated! 7/3/05
    enddo
  enddo
end SUBROUTINE contma

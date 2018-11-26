SUBROUTINE qdct1(shl,e,pois,c)
  use globalvar 
  implicit none
  include 'mpif.h'
  !
  !### program to read material property.
  !	B.D. 4/19/09
  !
  !...element arrays
  real (kind=8),dimension(nrowsh,nen) :: shl
   !...material properties
  real (kind=8),dimension(numat):: e,pois,rho,rdampm,rdampk
  real (kind=8),dimension(nrowc,nrowc,numat) :: c
  !...MPI
  integer (kind=4) :: master,me
  !
  !nee = nen * ned  !already in globalvar given explicitly
  !...compute bilinear shape functions and their local derivatives
  !		at integration points
  call qdcshl(shl)
   !...read material property data for 3D elastic continuum
  !		and compute material mouli matrix
  call prop3d(e,pois,c)
  !
end SUBROUTINE qdct1

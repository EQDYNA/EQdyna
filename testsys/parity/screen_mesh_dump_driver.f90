! Stand-alone driver for src/fortran/screenMesh.f90 (board row 167, step 1).
! Reads ONE mesh from a raw stream dump written by
! testsys/regression/test_mesh_screening.py (the planted-defect tests) and calls
! the production screenMesh on it, linked against the production globalvar.o,
! errorCodes.o and screenMesh.o. There is no hook in the production binary:
! the planted defect lives in the dump, and the verdict is this process's exit
! status (abortRun -> exit(code); 0 on a pass).
!
! Dump layout (native little-endian, in order):
!   nnode i8, nel i8, nx i4, ny i4, nz i4, nft i8, insertFaultType i4, tol f8
!   xline(nx) yline(ny) zline(nz) f8
!   meshCoor(3,nnode) f8, nodeElemIdRelation(8,nel) i8, elemTypeArr(nel) i4
!   nsmp(2,nft) i8, un(3,nft) us(3,nft) ud(3,nft) f8, arn(nft) f8
program screenMeshDumpDriver
    use globalvar
    implicit none
    character (len = 1024) :: path
    integer (kind = 8) :: nnode, nel, nft
    integer (kind = 4) :: nx, ny, nz, ift
    integer (kind = 4), allocatable :: nftnd0(:)
    real (kind = dp), allocatable :: xline(:), yline(:), zline(:)
    integer :: u

    if (command_argument_count() /= 1) then
        write(*,*) 'usage: screen_mesh_dump_driver <dump>'
        call exit(2)
    endif
    call get_command_argument(1, path)
    open(newunit=u, file=trim(path), access='stream', form='unformatted', status='old', action='read')
    read(u) nnode, nel, nx, ny, nz, nft, insertFaultType, tol
    ntotft = 1
    meshOnly = .true.
    allocate(xline(nx), yline(ny), zline(nz), nftnd0(1))
    nftnd0(1) = int(nft, 4)
    read(u) xline, yline, zline
    allocate(meshCoor(3,nnode), nodeElemIdRelation(8,nel), elemTypeArr(nel))
    read(u) meshCoor, nodeElemIdRelation, elemTypeArr
    allocate(nsmp(2,nft,1), un(3,nft,1), us(3,nft,1), ud(3,nft,1), arn(nft,1))
    read(u) nsmp, un, us, ud, arn
    close(u)
    call screenMesh(nnode, nel, nx, ny, nz, xline, yline, zline, nftnd0)
    write(*,'(a)') 'driver: screenMesh returned (no abort)'
end program screenMeshDumpDriver

! Stand-alone driver for src/fortran/screenMesh.f90 (board row 167, step 1).
! Reads ONE mesh from a raw stream dump written by
! testsys/regression/test_mesh_screening.py and calls the production screenMesh,
! linked against the production globalvar.o, errorCodes.o and screenMesh.o (no
! hook in the production binary). The verdict is this process's exit status
! (abortRun -> exit(code); 0 on a pass).
! Dump layout (native little-endian): nnode i8, nel i8, nft i8, tol f8, dx f8,
!   meshCoor(3,nnode) f8, nodeElemIdRelation(8,nel) i8, elemTypeArr(nel) i4,
!   nsmp(2,nft) i8, un(3,nft) f8
program screenMeshDumpDriver
    use globalvar
    implicit none
    character (len = 1024) :: path
    integer (kind = 8) :: nnode, nel, nft
    integer (kind = 4), allocatable :: nftnd0(:)
    integer :: u

    if (command_argument_count() /= 1) then
        write(*,*) 'usage: screen_mesh_dump_driver <dump>'
        call exit(2)
    endif
    call get_command_argument(1, path)
    open(newunit=u, file=trim(path), access='stream', form='unformatted', status='old', action='read')
    read(u) nnode, nel, nft, tol, dx
    ntotft = 1
    meshOnly = .true.
    allocate(nftnd0(1), meshCoor(3,nnode), nodeElemIdRelation(8,nel), elemTypeArr(nel), nsmp(2,nft,1), un(3,nft,1))
    nftnd0(1) = int(nft, 4)
    read(u) meshCoor, nodeElemIdRelation, elemTypeArr, nsmp, un
    close(u)
    call screenMesh(nnode, nel, nftnd0)
    write(*,'(a)') 'driver: screenMesh returned (no abort)'
end program screenMeshDumpDriver

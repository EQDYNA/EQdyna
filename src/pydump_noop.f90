! No-op stand-in for pydump.f90, linked into the DEFAULT (production) build.
! Keeps eqdyna3d.o/driver.o's compiled bytes identical between the default
! and PYDUMP=1 builds -- see driver.f90's comment on this hook.
subroutine pydump_state
    implicit none
end subroutine pydump_state

subroutine pydump_step(step)
    implicit none
    integer (kind = 4) :: step
end subroutine pydump_step

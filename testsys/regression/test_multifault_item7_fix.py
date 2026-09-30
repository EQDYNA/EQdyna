#! /usr/bin/env python3
"""
Regression guard (PROJECT_RULES.md rule 10; pathway_forward.md item 7) for
`netcdf_read_on_fault_eqdyna` writing every fault's on-fault netcdf input into
fault 1's `fric(...)` slots instead of its own fault's slots.

Background: `src/fortran/netcdf_io.f90:70-106` loops `do ift = 1, ntotft` /
`do i = 1, nftnd(ift)`, correctly using `ift` to index into `fxmin/fzmin` (the
per-fault grid origin used to look up `on_fault_vars(ii,jj,:)`) -- but every
`fric(FRIC_SLOT_*, i, 1)` assignment in that block used a LITERAL `1` for the
fault-index (3rd) dimension instead of `ift`. So on a hypothetical `ntotft>=2`
run, fault 2's on-fault input would overwrite fault 1's `fric` entries (fault
1's own data lost) and fault 2's own `fric` slots would be left at whatever
they were initialized to (never written). The restart routine right below it
(`netcdf_io.f90:159-179`, `netcdf_read_on_fault_eqdyna_restart`) already used
`ift` correctly throughout and was the template for this fix.

UPDATE (Row 17, multi-fault landing): `scripts/case.setup`'s ntotft>1 refusal
is gone (see test_multifault_two_fault_smoke.py's end-to-end 2-fault case),
and netcdf_io.f90 was generalized further in the same change: rather than one
shared (dip,strike) on-fault-variable grid read for every fault via per-fault
(ii,jj) offsets, each fault now gets its OWN full variable set in the netCDF
file, named with faultTag() (`''` for fault 1, `'ft2_'` for fault 2, ...) and
sized to that fault's own box. This test's synthetic netCDF file and expected
values were updated to match (two distinct variable sets, not one shared
one) -- the ORIGINAL item7 bug (fric indexed by literal `1` instead of `ift`)
is still exactly what is pinned; it is just no longer latent.

This test isolates the one subroutine directly: it compiles a tiny
standalone driver against the REAL, unmodified `netcdf_io.f90` (plus
`globalvar.f90`/`errorCodes.f90`/`library_output.f90`, its dependencies --
the last for the external `faultTag` function) with `mpif90` (the routine's
`use errorCodes` pulls in `abortRun`, which uses `include 'mpif.h'`), feeds
it a synthetic 2-fault `on_fault_vars_input.nc` (built with the real,
unmodified read-path, not a re-implementation) and 2 synthetic fault-node
placements, and asserts each fault's `fric` slots hold ITS OWN data, not the
other fault's (and not a third fault's-worth of aliasing onto the wrong
netCDF variable set either).

Fault-node placement is deliberately DIAGONAL in the (ii,jj) netCDF-variable
grid (fault 1's node -> cell (1,1); fault 2's node -> cell (2,2)) so the test
does not depend on getting netCDF-Fortran's dimension-reversal convention
(documented at the top of netcdf_io.f90) exactly right when building the
synthetic file from Python -- a diagonal cell's value is identical under
transpose, so any transpose ambiguity is moot here.

Cheap (rule 9): compiles 3 small .f90 files with mpif90 (no case, no mesh,
no actual MPI communication -- MPI_Initialized/MPI_Finalized inside abortRun
are never reached because the success path never calls abortRun) and runs a
single tiny netCDF read; well under a second. Fails loudly (rule 2) if
mpif90/netCDF-Fortran are unavailable -- never silently skips.
"""
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FSRC = os.path.join(ROOT, 'src', 'fortran')

sys.path.insert(0, ROOT)
from testsys.common import make_var  # noqa: E402

FORTRAN_DEPS = ('globalvar.f90', 'errorCodes.f90', 'library_output.f90', 'netcdf_io.f90')
# library_output.f90 added (row 17, multi-fault): netcdf_read_on_fault_eqdyna
# now calls the external function faultTag(ift), defined there.

NETCDF_INC = make_var('NETCDF_INC')
NETCDF_LIB = make_var('NETCDF_LIB')

# The 24 on-fault variable names netcdf_read_on_fault_eqdyna reads, in the
# EXACT order the subroutine's own nf90_inq_varid calls list them (this is
# what fixes each variable's slot index k = 1..24, used below to compute the
# expected fric value at each diagonal cell).
VAR_NAMES = [
    'sw_fs', 'sw_fd', 'sw_D0', 'rsf_a', 'rsf_b', 'rsf_Dc', 'rsf_v0', 'rsf_r0',
    'rsf_fw', 'rsf_vw', 'tp_a_hy', 'tp_a_th', 'tp_rouc', 'tp_lambda', 'tp_h',
    'tp_Tini', 'tp_pini', 'init_slip_rate', 'init_strike_shear',
    'init_normal_stress', 'init_state', 'tw_t0', 'cohesion', 'init_dip_shear',
]

# Driver: sets up two synthetic faults (ntotft=2), one fault-node pair each,
# placed so fault 1's node reads netCDF cell (1,1) and fault 2's node reads
# cell (2,2), then calls the REAL netcdf_read_on_fault_eqdyna and prints the
# resulting fric(FRIC_SLOT_SW_FS,1,:) and fric(FRIC_SLOT_COHESION,1,:) for
# both faults.
_DRIVER_SRC = r'''
program netcdf_item7_driver
    use globalvar
    implicit none

    ntotft = 2
    dx = 1.0d0
    dz = 1.0d0

    allocate(fxmin(2), fxmax(2), fzmin(2), fzmax(2))
    fxmin(1) = 0.0d0; fzmin(1) = 0.0d0; fxmax(1) = 1.0d0; fzmax(1) = 1.0d0
    fxmin(2) = 0.0d0; fzmin(2) = 0.0d0; fxmax(2) = 1.0d0; fzmax(2) = 1.0d0

    allocate(nftnd(2))
    nftnd(1) = 1
    nftnd(2) = 1

    allocate(nsmp(2,1,2))
    nsmp(1,1,1) = 1   ! fault 1's node pair -> mesh node id 1
    nsmp(1,1,2) = 2   ! fault 2's node pair -> mesh node id 2

    allocate(meshCoor(3,2))
    meshCoor(1,1) = 0.0d0; meshCoor(3,1) = 0.0d0   ! -> netCDF cell (1,1)
    meshCoor(1,2) = 1.0d0; meshCoor(3,2) = 1.0d0   ! -> netCDF cell (2,2)

    allocate(fric(100,1,2))
    fric = -999.0d0   ! sentinel: any slot the routine never touches stays -999

    call netcdf_read_on_fault_eqdyna

    write(*,'(A,ES25.16E3)') 'SW_FS_FAULT1 ', fric(FRIC_SLOT_SW_FS, 1, 1)
    write(*,'(A,ES25.16E3)') 'SW_FS_FAULT2 ', fric(FRIC_SLOT_SW_FS, 1, 2)
    write(*,'(A,ES25.16E3)') 'COHESION_FAULT1 ', fric(FRIC_SLOT_COHESION, 1, 1)
    write(*,'(A,ES25.16E3)') 'COHESION_FAULT2 ', fric(FRIC_SLOT_COHESION, 1, 2)
end program netcdf_item7_driver
'''


def build_netcdf_input(path):
    """Write a synthetic on_fault_vars_input.nc with 24 UNPREFIXED variables
    (fault 1's, per faultTag()'s convention: fault 1 takes no tag) and 24
    'ft2_'-prefixed variables (fault 2's) -- the multi-fault netCDF layout
    netcdf_write_on_fault_vars (scripts/case.setup) and
    netcdf_read_on_fault_eqdyna (src/fortran/netcdf_io.f90) now both use: one
    full variable set PER FAULT, named with faultTag(). Each is a 2x2 grid,
    diagonal cell (1,1)[0-based (0,0)] = 100+k and diagonal cell (2,2)
    [0-based (1,1)] = 200+k for fault 1's set, 300+k/400+k for fault 2's set
    (DIFFERENT numbers than fault 1's, so a regression that reads fault 1's
    variable for fault 2 -- old single-shared-grid behaviour -- is caught,
    not just the original item7 bug of writing fault 2's read into fault 1's
    fric slot). Off-diagonal cells are -1 (never read by this test's
    fault-node placement)."""
    import netCDF4
    import numpy as np

    with netCDF4.Dataset(path, 'w', format='NETCDF4') as ds:
        ds.createDimension('z', 2)
        ds.createDimension('x', 2)
        ds.createDimension('ft2_z', 2)
        ds.createDimension('ft2_x', 2)
        for k, name in enumerate(VAR_NAMES, start=1):
            var1 = ds.createVariable(name, 'f8', ('z', 'x'))
            data1 = np.full((2, 2), -1.0, dtype=np.float64)
            data1[0, 0] = 100.0 + k
            data1[1, 1] = 200.0 + k
            var1[:, :] = data1

            var2 = ds.createVariable('ft2_' + name, 'f8', ('ft2_z', 'ft2_x'))
            data2 = np.full((2, 2), -1.0, dtype=np.float64)
            data2[0, 0] = 300.0 + k
            data2[1, 1] = 400.0 + k
            var2[:, :] = data2


def build_driver(tmp):
    # ifort ignores gfortran's -J and writes .mod into the cwd (the repo root); cwd=tmp keeps them here.
    objs = []
    for fname in FORTRAN_DEPS:
        src = os.path.join(FSRC, fname)
        obj = os.path.join(tmp, fname.replace('.f90', '.o'))
        cmd = ['mpif90', '-O0', '-ffree-line-length-none', '-I', tmp, '-J', tmp]
        if fname == 'netcdf_io.f90':
            cmd += NETCDF_INC
        cmd += ['-c', src, '-o', obj]
        r = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True, timeout=60)
        if r.returncode != 0:
            raise RuntimeError(f'compiling {fname} failed:\n{r.stdout}\n{r.stderr}')
        objs.append(obj)

    driver_src = os.path.join(tmp, 'netcdf_item7_driver.f90')
    with open(driver_src, 'w') as fh:
        fh.write(_DRIVER_SRC)
    driver_obj = os.path.join(tmp, 'netcdf_item7_driver.o')
    r = subprocess.run(
        ['mpif90', '-O0', '-ffree-line-length-none', '-I', tmp, '-J', tmp,
         '-c', driver_src, '-o', driver_obj],
        cwd=tmp, capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        raise RuntimeError(f'compiling driver failed:\n{r.stdout}\n{r.stderr}')

    binary = os.path.join(tmp, 'netcdf_item7_driver')
    r = subprocess.run(
        ['mpif90', '-O0'] + objs + [driver_obj, '-o', binary] + NETCDF_LIB,
        cwd=tmp, capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        raise RuntimeError(f'linking driver failed:\n{r.stdout}\n{r.stderr}')
    return binary


def parse_output(text):
    vals = {}
    for line in text.strip().splitlines():
        parts = line.split()
        if len(parts) >= 2:
            vals[parts[0]] = float(parts[-1])
    return vals


def main():
    if shutil.which('mpif90') is None:
        print('FAIL test_multifault_item7_fix: mpif90 not on PATH '
              '(this guard requires a real Fortran compile, it does not skip silently)')
        return 1

    tmp = tempfile.mkdtemp(prefix='testMultifaultItem7.')
    try:
        try:
            binary = build_driver(tmp)
        except RuntimeError as e:
            print('FAIL test_multifault_item7_fix')
            print(' -', e)
            return 1

        build_netcdf_input(os.path.join(tmp, 'on_fault_vars_input.nc'))

        r = subprocess.run([binary], cwd=tmp, capture_output=True, text=True, timeout=30)
        if r.returncode != 0:
            print('FAIL test_multifault_item7_fix')
            print(f' - driver exited {r.returncode}:\n{r.stdout}\n{r.stderr}')
            return 1
        vals = parse_output(r.stdout)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    expected = {
        'SW_FS_FAULT1': 101.0,     # 100 + slot 1 (sw_fs), unprefixed set
        'SW_FS_FAULT2': 401.0,     # 400 + slot 1 (sw_fs), 'ft2_' set
        'COHESION_FAULT1': 123.0,  # 100 + slot 23 (cohesion), unprefixed set
        'COHESION_FAULT2': 423.0,  # 400 + slot 23 (cohesion), 'ft2_' set
    }
    fails = []
    for key, exp in expected.items():
        got = vals.get(key)
        if got is None:
            fails.append(f'{key}: missing from driver output ({vals!r})')
        elif abs(got - exp) > 1e-6:
            fails.append(f'{key}: expected {exp!r}, got {got!r} '
                          '(fault-2 data overwriting fault-1 slot, fault-2 slot never written, '
                          'or fault 2 reading the UNPREFIXED set instead of its own ft2_ set)')

    # Distinctness check: if the original item7 bug were back, both faults
    # would collapse onto the SAME slot (fault 1's), so fault2's own numbers
    # would equal -999 (never written) rather than 401/423. If instead fault
    # 2 fell back to reading the unprefixed ('sw_fs') set -- the
    # single-shared-grid behaviour this test used to assume -- its values
    # would equal fault 1's (101/123) rather than 401/423.
    if vals.get('SW_FS_FAULT2') == -999.0 or vals.get('COHESION_FAULT2') == -999.0:
        fails.append('fault 2 fric slots were never written (bug: literal `1` '
                      'still used instead of `ift`)')
    if vals.get('SW_FS_FAULT2') == vals.get('SW_FS_FAULT1'):
        fails.append('fault 2 read the UNPREFIXED (fault 1) netCDF variable set '
                      'instead of its own ft2_-prefixed set')

    if fails:
        print('FAIL test_multifault_item7_fix')
        for f in fails:
            print('  -', f)
        return 1
    print('SUCCESS test_multifault_item7_fix '
          '(fault 1 and fault 2 fric slots independently populated from their own on-fault netCDF data)')
    return 0


if __name__ == '__main__':
    sys.exit(main())

#! /usr/bin/env python3
"""
Regression guard (PROJECT_RULES.md rule 10; pathway_forward.md item 10,
superseded by Row 17's multi-fault landing) for `MPI4arn`
(`src/fortran/meshgen.f90`) corrupting one fault's per-boundary state
(`fltgm`/`fltl`/`fltr`/`fltf`/`fltb`/`fltd`/`fltu`/`fltnum`) when a LATER
fault is processed.

HISTORY: item 10's original defect was `MPI4arn` allocating `fltl..fltu` from
the PREVIOUS fault's leftover `fltnum` counts without deallocating an array a
previous call already allocated -- a bug in a single-fault-sized,
allocate-on-every-call design. Row 17 replaced that design entirely: these
arrays are now `(nftmx, ntotft)` / `(6, ntotft)`, pre-allocated ONCE
(`eqdyna3d.f90`'s `allocInit`) and filled in place, one COLUMN per fault, by
`MPI4arn` (no allocate/deallocate inside it at all any more -- see that
subroutine's own header comment for why). The ORIGINAL failure mode (fault 2
corrupting fault 1's already-recorded state) is exactly what is still worth
guarding: with a shared column layout, a bug that indexes the wrong column,
or writes past a column's own bounds, would have the same symptom (fault 1's
state changes when fault 2 is processed) even though the mechanism (stale
`allocate`) is gone.

This test isolates `MPI4arn` directly against the REAL, unmodified
`meshgen.f90` (and everything it links against -- built exactly as
production, via the project's own `make`, no dependency hand-picked or
reimplemented), calling it twice in sequence to simulate a 2-fault
`do ift=1,ntotft` loop:
  - fault 1: 3 local fault-node-pairs on boundaries left/right/down
  - fault 2: 2 local fault-node-pairs, BOTH on the left boundary -- a
    DIFFERENT count than fault 1's (1 -> 2), so a correct implementation must
    both leave fault 1's column alone AND size fault 2's own column
    correctly, not just happen to not crash.
npx=npy=npz=1 so the MPI send/recv branches (which need a real multi-rank
run) are never entered -- this isolates exactly the per-fault bookkeeping,
not the (already separately guarded) cross-rank arn exchange.

The genuine end-to-end case (test_multifault_two_fault_smoke.py) additionally
exercises the cross-rank exchange for two real faults; this test is cheaper
and pins the bookkeeping in isolation.

Cheap-ish (rule 9): one `make` of the 26-file src/fortran tree (already-built
objects are reused if the tree hasn't changed) plus one tiny compile/link/run
of a driver that does no meshing, no case, no actual MPI communication.
Fails loudly (rule 2) if mpif90/the src/ build are unavailable -- never
silently skips.
"""
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FSRC = os.path.join(ROOT, 'src', 'fortran')

sys.path.insert(0, ROOT)
from testsys.common import MACHINE, make_var  # noqa: E402

# Every production object file EXCEPT eqdyna3d.o (which carries `program
# EQdyna` -- our driver supplies its own program unit instead). Confirmed by
# inspection that no other file calls any subroutine defined only in
# eqdyna3d.f90 (allocInit, allocInitAfterMeshGen, init_vel,
# checkFaultMPIAlignment, checkMeshMaterial, checkArrSize), so dropping
# eqdyna3d.o cannot leave an unresolved symbol.
OBJS = [
    'globalvar.o', 'errorCodes.o', 'countMeshEntities.o', 'meshgen.o',
    'calcLocalShapeFunc.o', 'driver.o', 'assembleGlobalMass.o',
    'calcGlobalShapeFunc.o', 'calcB.o', 'assembleGlobalKU.o',
    'calcElemMass.o', 'calcElemKU.o', 'calcHourglassResist.o', 'fric.o',
    'faulting.o', 'computePMLDampingVector.o', 'calcQAttenuationCoeff.o',
    'readInputFiles.o', 'updateThermalPressurization.o', 'library.o',
    'func_lib.o', 'checkInputConsistency.o', 'library_degeneration.o',
    'library_output.o', 'netcdf_io.o',
]

NETCDF_LIB = make_var('NETCDF_LIB')

# Driver: simulates meshgen's `do ift=1,ntotft: call MPI4arn(...)` loop for
# two faults back to back, replicating the (nftmx, ntotft) pre-allocation
# eqdyna3d.f90's allocInit does in production (MPI4arn no longer allocates
# fltl..fltu itself -- it fills a pre-sized column in place).
_DRIVER_SRC = r'''
program mpi4arn_item10_driver
    use globalvar
    implicit none
    integer, parameter :: NLOCAL = 10

    ntotft = 2
    npx = 1
    npy = 1
    npz = 1
    me = 0

    allocate(fltgm(NLOCAL,ntotft), fltl(NLOCAL,ntotft), fltr(NLOCAL,ntotft), &
             fltf(NLOCAL,ntotft), fltb(NLOCAL,ntotft), fltd(NLOCAL,ntotft), &
             fltu(NLOCAL,ntotft), fltnum(6,ntotft))
    fltgm = 0
    fltl = 0; fltr = 0; fltf = 0; fltb = 0; fltd = 0; fltu = 0
    fltnum = 0
    fltMPI = .false.

    ! ---- Fault 1: 3 local fault-node-pairs -- left(1), right(2), down(100).
    fltgm(1,1) = 1
    fltgm(2,1) = 2
    fltgm(3,1) = 100

    call MPI4arn(1, 1, 1, 0, 0, 0, 3, 1)

    print *, 'FLTNUM1', fltnum(:,1)
    print *, 'FLTL1', fltl(:,1)
    print *, 'FLTR1', fltr(:,1)
    print *, 'FLTD1', fltd(:,1)

    ! ---- Fault 2: 2 local fault-node-pairs, BOTH on the left boundary (a
    ! DIFFERENT fltnum(1,2) than fault 1's fltnum(1,1) -- 1 -> 2 -- to prove
    ! correct per-fault sizing, not just reuse of a coincidentally-sized
    ! shared array).
    fltgm(1,2) = 1
    fltgm(2,2) = 1

    call MPI4arn(1, 1, 1, 0, 0, 0, 2, 2)

    print *, 'FLTNUM1_AFTER', fltnum(:,1)
    print *, 'FLTL1_AFTER', fltl(:,1)
    print *, 'FLTR1_AFTER', fltr(:,1)
    print *, 'FLTD1_AFTER', fltd(:,1)
    print *, 'FLTNUM2', fltnum(:,2)
    print *, 'FLTL2', fltl(:,2)
    print *, 'DRIVER_OK'
end program mpi4arn_item10_driver
'''


def build_production_objects():
    env = dict(os.environ)
    env['MACHINE'] = MACHINE
    r = subprocess.run(['make'], cwd=FSRC, env=env, capture_output=True, text=True, timeout=300)
    missing = [o for o in OBJS if not os.path.exists(os.path.join(FSRC, o))]
    if r.returncode != 0 or missing:
        raise RuntimeError(
            f'src/fortran build failed (exit {r.returncode}) or missing objects {missing}; '
            f'is mpif90 (MACHINE={MACHINE}) on PATH?\n' + r.stdout[-3000:] + r.stderr[-3000:])


def build_driver(tmp):
    driver_src = os.path.join(tmp, 'mpi4arn_item10_driver.f90')
    with open(driver_src, 'w') as fh:
        fh.write(_DRIVER_SRC)
    driver_obj = os.path.join(tmp, 'mpi4arn_item10_driver.o')
    r = subprocess.run(
        ['mpif90', '-O0', '-ffree-line-length-none', '-I', FSRC, '-J', tmp,
         '-c', driver_src, '-o', driver_obj],
        capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        raise RuntimeError(f'compiling driver failed:\n{r.stdout}\n{r.stderr}')

    binary = os.path.join(tmp, 'mpi4arn_item10_driver')
    objs = [os.path.join(FSRC, o) for o in OBJS]
    r = subprocess.run(
        ['mpif90', '-fopenmp', '-ffree-line-length-none', '-O0'] + objs +
        [driver_obj, '-o', binary] + NETCDF_LIB,
        capture_output=True, text=True, timeout=120)
    if r.returncode != 0:
        raise RuntimeError(f'linking driver failed:\n{r.stdout}\n{r.stderr}')
    return binary


def parse(text):
    """Return {marker: [ints]} for each 'print *, MARKER v1 v2 ...' line."""
    out = {}
    for line in text.splitlines():
        parts = line.split()
        if not parts:
            continue
        marker = parts[0]
        rest = parts[1:]
        out[marker] = rest
    return out


def as_ints(vals):
    return [int(v) for v in vals]


def main():
    if shutil.which('mpif90') is None:
        print('FAIL test_multifault_item10_fix: mpif90 not on PATH '
              '(this guard requires a real Fortran build, it does not skip silently)')
        return 1

    try:
        build_production_objects()
    except RuntimeError as e:
        print('FAIL test_multifault_item10_fix')
        print(' -', e)
        return 1

    tmp = tempfile.mkdtemp(prefix='testMultifaultItem10.')
    try:
        try:
            binary = build_driver(tmp)
        except RuntimeError as e:
            print('FAIL test_multifault_item10_fix')
            print(' -', e)
            return 1

        r = subprocess.run([binary], cwd=tmp, capture_output=True, text=True, timeout=30)
        if r.returncode != 0:
            print('FAIL test_multifault_item10_fix')
            print(f' - driver exited {r.returncode} (fault-2 call likely aborted):\n{r.stdout}\n{r.stderr}')
            return 1
        vals = parse(r.stdout)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    fails = []
    if 'DRIVER_OK' not in vals:
        fails.append(f'driver did not report DRIVER_OK -- output:\n{r.stdout}')

    # After fault 1's call: fltnum(1)=1 (left), fltnum(2)=1 (right),
    # fltnum(5)=1 (down); fltl/fltr/fltd each hold local node index 1/2/3.
    want_fltnum1 = [1, 1, 0, 0, 1, 0]
    fltnum1 = as_ints(vals.get('FLTNUM1', []))
    if fltnum1 != want_fltnum1:
        fails.append(f'FLTNUM1: expected {want_fltnum1}, got {fltnum1}')
    if as_ints(vals.get('FLTL1', [])) != [1,0,0,0,0,0,0,0,0,0]:
        fails.append(f'FLTL1: expected [1,0,...], got {vals.get("FLTL1")}')
    if as_ints(vals.get('FLTR1', [])) != [2,0,0,0,0,0,0,0,0,0]:
        fails.append(f'FLTR1: expected [2,0,...], got {vals.get("FLTR1")}')
    if as_ints(vals.get('FLTD1', [])) != [3,0,0,0,0,0,0,0,0,0]:
        fails.append(f'FLTD1: expected [3,0,...], got {vals.get("FLTD1")}')

    # THE ACTUAL GUARD: after fault 2's call, fault 1's OWN column must be
    # UNCHANGED -- this is what a cross-column corruption (wrong index into
    # the shared (nftmx, ntotft) array, or an out-of-bounds write from
    # fault 2's fill loop) would break.
    fltnum1_after = as_ints(vals.get('FLTNUM1_AFTER', []))
    if fltnum1_after != want_fltnum1:
        fails.append(f'FLTNUM1_AFTER (fault 1''s column, AFTER fault 2''s call): '
                      f'expected {want_fltnum1} (unchanged), got {fltnum1_after} '
                      '-- fault 2 corrupted fault 1''s state')
    if as_ints(vals.get('FLTL1_AFTER', [])) != [1,0,0,0,0,0,0,0,0,0]:
        fails.append(f'FLTL1_AFTER: expected [1,0,...] (unchanged), got {vals.get("FLTL1_AFTER")}')
    if as_ints(vals.get('FLTR1_AFTER', [])) != [2,0,0,0,0,0,0,0,0,0]:
        fails.append(f'FLTR1_AFTER: expected [2,0,...] (unchanged), got {vals.get("FLTR1_AFTER")}')
    if as_ints(vals.get('FLTD1_AFTER', [])) != [3,0,0,0,0,0,0,0,0,0]:
        fails.append(f'FLTD1_AFTER: expected [3,0,...] (unchanged), got {vals.get("FLTD1_AFTER")}')

    # Fault 2's own column: its OWN count (2, not fault 1's leftover 1), and
    # its own local node indices (1, 2).
    want_fltnum2 = [2, 0, 0, 0, 0, 0]
    fltnum2 = as_ints(vals.get('FLTNUM2', []))
    if fltnum2 != want_fltnum2:
        fails.append(f'FLTNUM2: expected {want_fltnum2} (fault 2 own counts), got {fltnum2}')
    if as_ints(vals.get('FLTL2', [])) != [1,2,0,0,0,0,0,0,0,0]:
        fails.append(f'FLTL2: expected [1,2,0,...] (both of fault 2\'s own nodes), got {vals.get("FLTL2")}')

    if fails:
        print('FAIL test_multifault_item10_fix')
        for f in fails:
            print('  -', f)
        return 1
    print('SUCCESS test_multifault_item10_fix '
          '(fault 2''s MPI4arn call left fault 1''s column untouched and filled its own column correctly)')
    return 0


if __name__ == '__main__':
    sys.exit(main())

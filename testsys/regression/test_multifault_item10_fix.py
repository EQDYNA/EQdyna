#! /usr/bin/env python3
"""
Regression guard (PROJECT_RULES.md rule 10; pathway_forward.md item 10) for
`MPI4arn` (`src/fortran/meshgen.f90`) allocating `fltl..fltu` from the
PREVIOUS fault's leftover `fltnum` counts, without deallocating an array a
previous call already allocated.

Background: `meshgen` calls `MPI4arn` once per fault, from `do ift=1,ntotft`.
On entry, `MPI4arn` used to do:

    if(fltnum(1) /= 0) allocate(fltl(fltnum(1)))   ! ... same for r/f/b/d/u
    fltnum = 0
    do i = 1, totalNumFaultNode                    ! recompute fltnum + fill
        ...
    enddo

`fltnum` on entry is stale: either createMasterNode's running tally (first
call) or -- for every call after the first -- exactly what THIS SAME loop
left it at for the PREVIOUS fault. So on fault 2's call, `fltnum(k)` from
fault 1 is still nonzero, `allocate(fltl(fltnum(1)))` fires again, and `fltl`
is already allocated from fault 1's own call: Fortran's `allocate` on an
already-allocated object is a runtime error. `ntotft==1` (the only case that
runs today) never hits this because there is only ever one call.

Fix: count THIS fault's own per-boundary membership from `fltgm` in a
pre-pass (the same data + arithmetic the fill loop below already reads/uses,
just counted before allocating instead of after), deallocate any array a
previous call left allocated, then (re)allocate to the freshly-counted size.

This bug is LATENT and cannot be exercised end-to-end (case.setup refuses
ntotft>1, pathway_forward.md item 17 -- deliberately left refused, not
touched here). This test isolates `MPI4arn` directly against the REAL,
unmodified `meshgen.f90` (and everything it links against -- built exactly
as production, via the project's own `make`, so no dependency is
hand-picked or reimplemented) and calls it twice in sequence, simulating a
2-fault `do ift=1,ntotft` loop:
  - fault 1: 3 local fault-node-pairs on boundaries left/right/down
  - fault 2: 2 local fault-node-pairs, BOTH on the left boundary -- a
    DIFFERENT count than fault 1's (1 -> 2), so a correct fix must both
    avoid the double-allocate AND resize to the new fault's own count,
    not just reuse-with-luck a same-sized stale array.
npx=npy=npz=1 so the MPI send/recv branches (which need a real multi-rank
run) are never entered -- this isolates exactly the allocate/reset ordering
bug, not the (already separately guarded) cross-rank arn exchange.

Verified against the pre-fix code directly (not shipped, sanity-checked by
hand during development of this test): the unfixed subroutine aborts on the
second call with a gfortran "Attempting to allocate already-allocated
variable" runtime error; the fixed subroutine returns fltnum/fltl/fltr/fltd
that reflect ONLY fault 2's own data after the second call.

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
MACHINE = os.environ.get('EQDYNA_TEST_MACHINE', 'ubuntu')

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

NETCDF_LIB = ['-L', '/usr/lib/x86_64-linux-gnu', '-lnetcdf', '-lnetcdff']

# Driver: simulates meshgen's `do ift=1,ntotft: call MPI4arn(...)` loop for
# two faults back to back. npx=npy=npz=1 so the cross-rank exchange branches
# (which need a real MPI run) are never entered.
_DRIVER_SRC = r'''
program mpi4arn_item10_driver
    use globalvar
    implicit none

    ntotft = 2
    npx = 1
    npy = 1
    npz = 1
    me = 0

    allocate(fltgm(10))
    fltgm = 0

    ! ---- Fault 1: 3 local fault-node-pairs -- left(1), right(2), down(100).
    fltgm(1) = 1
    fltgm(2) = 2
    fltgm(3) = 100
    ! Pre-set fltnum as createMasterNode's running tally would have left it
    ! by the time meshgen's own MPI4arn(ift=1) call happens in the real flow.
    fltnum = 0
    fltnum(1) = 1
    fltnum(2) = 1
    fltnum(5) = 1

    call MPI4arn(1, 1, 1, 0, 0, 0, 3, 1)

    print *, 'FLTNUM1', fltnum
    print *, 'FLTL1', fltl
    print *, 'FLTR1', fltr
    print *, 'FLTD1', fltd

    ! ---- Fault 2: 2 local fault-node-pairs, BOTH on the left boundary (a
    ! DIFFERENT fltnum(1) than fault 1's -- 1 -> 2 -- to prove correct
    ! resizing, not just reuse of a stale-but-coincidentally-sized array).
    fltgm(1) = 1
    fltgm(2) = 1

    call MPI4arn(1, 1, 1, 0, 0, 0, 2, 2)

    print *, 'FLTNUM2', fltnum
    print *, 'FLTL2', fltl
    print *, 'FLTR2_ALLOCATED', allocated(fltr)
    print *, 'FLTD2_ALLOCATED', allocated(fltd)
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
            print(f' - driver exited {r.returncode} (fault-2 call likely aborted -- '
                  f'the double-allocate bug is back):\n{r.stdout}\n{r.stderr}')
            return 1
        vals = parse(r.stdout)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    fails = []
    if 'DRIVER_OK' not in vals and not any(k == 'DRIVER_OK' for k in vals):
        fails.append(f'driver did not report DRIVER_OK -- output:\n{r.stdout}')

    # After fault 1's call: fltnum(1)=1 (left), fltnum(2)=1 (right),
    # fltnum(5)=1 (down); fltl/fltr/fltd each hold local node index 1/2/3.
    fltnum1 = as_ints(vals.get('FLTNUM1', []))
    if fltnum1 != [1, 1, 0, 0, 1, 0]:
        fails.append(f'FLTNUM1: expected [1, 1, 0, 0, 1, 0], got {fltnum1}')
    if as_ints(vals.get('FLTL1', [])) != [1]:
        fails.append(f'FLTL1: expected [1], got {vals.get("FLTL1")}')
    if as_ints(vals.get('FLTR1', [])) != [2]:
        fails.append(f'FLTR1: expected [2], got {vals.get("FLTR1")}')
    if as_ints(vals.get('FLTD1', [])) != [3]:
        fails.append(f'FLTD1: expected [3], got {vals.get("FLTD1")}')

    # After fault 2's call: the whole point of the fix. fltnum(1) must be 2
    # (fault 2's OWN count, not fault 1's leftover 1), and fltl must hold
    # fault 2's own local node indices (1, 2), not fault 1's stale [1].
    fltnum2 = as_ints(vals.get('FLTNUM2', []))
    if fltnum2 != [2, 0, 0, 0, 0, 0]:
        fails.append(f'FLTNUM2: expected [2, 0, 0, 0, 0, 0] (fault 2 own counts), got {fltnum2} '
                      '(stale fault-1 counts would look like [1, 1, 0, 0, 1, 0])')
    if as_ints(vals.get('FLTL2', [])) != [1, 2]:
        fails.append(f'FLTL2: expected [1, 2] (both of fault 2\'s own nodes), got {vals.get("FLTL2")}')
    # fault 2 has zero right/down-boundary nodes -- a correct fix deallocates
    # fault 1's leftover fltr/fltd rather than leaving them allocated (and
    # sized for fault 1) with a now-zero fltnum guarding their use.
    if vals.get('FLTR2_ALLOCATED') != ['F']:
        fails.append(f'FLTR2_ALLOCATED: expected F (deallocated, fault 2 has none), '
                      f'got {vals.get("FLTR2_ALLOCATED")}')
    if vals.get('FLTD2_ALLOCATED') != ['F']:
        fails.append(f'FLTD2_ALLOCATED: expected F (deallocated, fault 2 has none), '
                      f'got {vals.get("FLTD2_ALLOCATED")}')

    if fails:
        print('FAIL test_multifault_item10_fix')
        for f in fails:
            print('  -', f)
        return 1
    print('SUCCESS test_multifault_item10_fix '
          '(fault 2 call did not crash and produced fault-2-only fltnum/fltl, not fault-1 leftovers)')
    return 0


if __name__ == '__main__':
    sys.exit(main())

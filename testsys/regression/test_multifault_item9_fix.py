#! /usr/bin/env python3
"""
Regression guard (PROJECT_RULES.md rule 10; pathway_forward.md item 9) for
`output_onfault_st` (`src/fortran/library_output.f90`) opening unit 51 only
for "main fault" (j==1) stations while writing to it unconditionally for
every station -- see item 9 for the incident this fixed.

The fix: the filename/header/open block now runs for EVERY station
regardless of `j = anonfs(3,i)` -- it is fully self-contained per loop
iteration (open, write, close(51) all within the same pass over station
`i`), so there was nothing "main-fault-only" about it. No-op for ntotft==1
(there `j` is always 1, so the block always ran before too).

This bug is LATENT and cannot be exercised end-to-end (case.setup refuses
ntotft>1, pathway_forward.md item 17 -- deliberately left refused, not
touched here). This test isolates the subroutine directly: `output_onfault_st`
depends on nothing but the `globalvar` module (no MPI, no other source
file), so it compiles clean with plain gfortran. The driver feeds it TWO
synthetic on-fault stations -- one on fault 1, one on fault 2 -- with
distinct strike/dip coordinates (so each gets its own, distinguishable,
correctly-named output file) and a distinct marker value in each station's
time series, then asserts BOTH expected files exist with the RIGHT station's
data in each.

Cheap (rule 9): compiles 2 small .f90 files with plain gfortran (no MPI, no
case, no mesh) and writes ~2 tiny text files; well under a second. Fails
loudly (rule 2) if gfortran is unavailable -- never silently skips.
"""
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FSRC = os.path.join(ROOT, 'src', 'fortran')

FORTRAN_DEPS = ('globalvar.f90', 'library_output.f90')

# Driver: two on-fault stations, station 1 on fault 1, station 2 on fault 2.
# Station 1 -> filename faultst005dp003.txt (strike 500 m -> 005, dip 300 m
# with fltxyz(2,4,1)=pi/2 -> dsin=1 -> 003).
# Station 2 -> filename faultst010dp004.txt (strike 1000 m -> 010, dip 400 m
# -> 004). Distinct from station 1's filename by construction.
# Each station's onFaultQuantHistSCECForm(1,1,i) (the Time column) carries a
# distinct marker (1.111 / 2.222) so the test can tell the two files apart
# by content, not just by existing.
_DRIVER_SRC = r'''
program library_output_item9_driver
    use globalvar
    implicit none

    ntotft = 2
    friclaw = 1
    nstep = 1
    dx = 100.0d0
    dt = 0.01d0
    numOfOnFaultStCount = 2

    allocate(anonfs(3,2))
    anonfs(1,1) = 1; anonfs(2,1) = 1; anonfs(3,1) = 1   ! station 1: fault 1
    anonfs(1,2) = 1; anonfs(2,2) = 1; anonfs(3,2) = 2   ! station 2: fault 2

    allocate(xonfs(2,1,2))
    xonfs(1,1,1) = 500.0d0;  xonfs(2,1,1) = -300.0d0    ! fault 1's station
    xonfs(1,1,2) = 1000.0d0; xonfs(2,1,2) = -400.0d0    ! fault 2's station

    allocate(fltxyz(2,4,2))
    fltxyz(2,4,1) = pi/2.0d0   ! dip = 90 deg, dsin = 1 (used for both, see item 9 note)
    fltxyz(2,4,2) = pi/2.0d0

    allocate(onFaultQuantHistSCECForm(12,1,2))
    onFaultQuantHistSCECForm = 0.0d0
    onFaultQuantHistSCECForm(1,1,1) = 1.111d0   ! station 1's marker
    onFaultQuantHistSCECForm(1,1,2) = 2.222d0   ! station 2's marker

    call output_onfault_st
end program library_output_item9_driver
'''


def build_driver(tmp):
    objs = []
    for fname in FORTRAN_DEPS:
        src = os.path.join(FSRC, fname)
        obj = os.path.join(tmp, fname.replace('.f90', '.o'))
        r = subprocess.run(
            ['gfortran', '-O0', '-ffree-line-length-none', '-I', tmp, '-J', tmp,
             '-c', src, '-o', obj],
            capture_output=True, text=True, timeout=60)
        if r.returncode != 0:
            raise RuntimeError(f'compiling {fname} failed:\n{r.stdout}\n{r.stderr}')
        objs.append(obj)

    driver_src = os.path.join(tmp, 'library_output_item9_driver.f90')
    with open(driver_src, 'w') as fh:
        fh.write(_DRIVER_SRC)
    driver_obj = os.path.join(tmp, 'library_output_item9_driver.o')
    r = subprocess.run(
        ['gfortran', '-O0', '-ffree-line-length-none', '-I', tmp, '-J', tmp,
         '-c', driver_src, '-o', driver_obj],
        capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        raise RuntimeError(f'compiling driver failed:\n{r.stdout}\n{r.stderr}')

    binary = os.path.join(tmp, 'library_output_item9_driver')
    r = subprocess.run(
        ['gfortran', '-O0'] + objs + [driver_obj, '-o', binary],
        capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        raise RuntimeError(f'linking driver failed:\n{r.stdout}\n{r.stderr}')
    return binary


def main():
    if shutil.which('gfortran') is None:
        print('FAIL test_multifault_item9_fix: gfortran not on PATH '
              '(this guard requires a real Fortran compile, it does not skip silently)')
        return 1

    tmp = tempfile.mkdtemp(prefix='testMultifaultItem9.')
    fails = []
    try:
        try:
            binary = build_driver(tmp)
        except RuntimeError as e:
            print('FAIL test_multifault_item9_fix')
            print(' -', e)
            return 1

        r = subprocess.run([binary], cwd=tmp, capture_output=True, text=True, timeout=30)
        if r.returncode != 0:
            print('FAIL test_multifault_item9_fix')
            print(f' - driver exited {r.returncode}:\n{r.stdout}\n{r.stderr}')
            return 1

        # The compiler-default fallback for an unconnected unit (gfortran:
        # fort.51) must never appear -- its presence means fault 2's station
        # skipped the open and fell through to writing an unconnected unit.
        stray = os.path.join(tmp, 'fort.51')
        if os.path.exists(stray):
            fails.append('fort.51 (gfortran default for an unopened unit) was created -- '
                          'fault 2 station wrote to an unconnected unit')

        expected = {
            os.path.join(tmp, 'faultst005dp003.txt'): '0.1111',   # fault 1's station
            os.path.join(tmp, 'faultst010dp004.txt'): '0.2222',   # fault 2's station
        }
        for path, marker in expected.items():
            if not os.path.exists(path):
                fails.append(f'{os.path.basename(path)} was never created '
                              f'(fault-2 station output missing or misdirected)')
                continue
            with open(path) as fh:
                content = fh.read()
            if marker not in content:
                fails.append(f'{os.path.basename(path)} does not contain its own marker '
                              f'{marker!r} -- got:\n{content}')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    if fails:
        print('FAIL test_multifault_item9_fix')
        for f in fails:
            print('  -', f)
        return 1
    print('SUCCESS test_multifault_item9_fix '
          '(fault 1 and fault 2 on-fault stations each got their own, correctly-named output file)')
    return 0


if __name__ == '__main__':
    sys.exit(main())

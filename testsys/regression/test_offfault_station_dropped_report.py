#! /usr/bin/env python3
"""
Regression guard (pathway_forward.md item 94) for off-fault stations that
match no grid node being dropped with NO message.

THE BUG: `setSurfaceStation` (`src/fortran/meshgen.f90`) matches a station's
depth EXACTLY (|nodeCoor(3) - x4nds(3,i)| < tol) and snaps only x and y to
the nearest interior node. A station whose depth is not a grid z-plane (or
that lies outside the interior x/y range) therefore matches no node on any
rank, and `output_offfault_st` writes no file for it -- silently. test.tpv8
at dx = 500 m requested 15 off-fault stations and wrote 11: its four
z = -0.3 km stations sit between the z = 0 and z = -0.5 km planes.

THE FIX: `checkOffFaultStationCoverage` (`eqdyna3d.f90`) OR-reduces `n4yn`
over all ranks after `meshgen`, and rank 0 calls
`report_dropped_offfault_st` (`library_output.f90`), which prints a WARNING
with the dropped count and one line per dropped station (index and x,y,z in
km). It neither snaps nor refuses; that choice changes which files a case
writes and is the owner's.

This test compiles `report_dropped_offfault_st` directly (no MPI, no case,
no mesh; same technique as test_station_header_column_count.py) and checks
BEHAVIOUR, not source text (rule 10a):
  (a) 3 stations, station 2 unmatched -> the WARNING says "1 of 3", names
      station 2 with its km coordinates, and names neither 1 nor 3;
  (b) all 3 matched -> prints nothing at all.
Before the fix the subroutine did not exist, so the build itself fails.
Fails loudly (rule 2) if gfortran is unavailable -- never silently skips.
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FSRC = os.path.join(ROOT, 'src', 'fortran')
FORTRAN_DEPS = ('globalvar.f90', 'library_output.f90')

_DRIVER_SRC = r'''
program item94_driver
    use globalvar
    implicit none
    logical :: matched(3)

    totalNumOfOffSt = 3
    allocate(x4nds(3,3))
    x4nds(:,1) = (/    0.0d0,  1000.0d0,    0.0d0 /)
    x4nds(:,2) = (/    0.0d0,  -500.0d0, -300.0d0 /)
    x4nds(:,3) = (/ 12000.0d0, 3000.0d0, -12000.0d0 /)
    matched = (/ {m1}, {m2}, {m3} /)
    call report_dropped_offfault_st(matched)
end program item94_driver
'''


def build_and_run(tmp, tag, flags):
    objs = []
    for fname in FORTRAN_DEPS:
        obj = os.path.join(tmp, fname.replace('.f90', f'.{tag}.o'))
        r = subprocess.run(
            ['gfortran', '-O0', '-ffree-line-length-none', '-I', tmp, '-J', tmp,
             '-c', os.path.join(FSRC, fname), '-o', obj],
            capture_output=True, text=True, timeout=60)
        if r.returncode != 0:
            raise RuntimeError(f'compiling {fname} ({tag}) failed:\n{r.stdout}\n{r.stderr}')
        objs.append(obj)
    src = os.path.join(tmp, f'driver.{tag}.f90')
    with open(src, 'w') as fh:
        fh.write(_DRIVER_SRC.format(m1=flags[0], m2=flags[1], m3=flags[2]))
    binary = os.path.join(tmp, f'driver.{tag}')
    r = subprocess.run(
        ['gfortran', '-O0', '-ffree-line-length-none', '-I', tmp, '-J', tmp]
        + objs + [src, '-o', binary],
        capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        raise RuntimeError(f'building driver ({tag}) failed:\n{r.stdout}\n{r.stderr}')
    r = subprocess.run([binary], cwd=tmp, capture_output=True, text=True, timeout=30)
    if r.returncode != 0:
        raise RuntimeError(f'driver ({tag}) exited {r.returncode}:\n{r.stdout}\n{r.stderr}')
    return r.stdout


def main():
    if shutil.which('gfortran') is None:
        print('FAIL test_offfault_station_dropped_report: gfortran not on PATH '
              '(this guard requires a real Fortran compile, it does not skip silently)')
        return 1
    tmp = tempfile.mkdtemp(prefix='testOffFaultDropped.')
    fails = []
    try:
        try:
            out = build_and_run(tmp, 'one', ('.true.', '.false.', '.true.'))
            quiet = build_and_run(tmp, 'none', ('.true.', '.true.', '.true.'))
        except RuntimeError as e:
            print('FAIL test_offfault_station_dropped_report')
            print(' -', e)
            return 1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    named = [int(n) for n in re.findall(r'dropped off-fault station (\d+)', out)]
    if not re.search(r'WARNING: 1 of 3 requested off-fault stations', out):
        fails.append(f'(a) no "WARNING: 1 of 3" line; stdout was:\n{out}')
    if named != [2]:
        fails.append(f'(a) named stations {named}, expected [2]; stdout was:\n{out}')
    m = re.search(r'station 2 at x,y,z =\s*(\S+)\s+(\S+)\s+(\S+)\s+km', out)
    if not m or [float(v) for v in m.groups()] != [0.0, -0.5, -0.3]:
        fails.append(f'(a) station 2 coordinates not reported as 0.000 -0.500 -0.300 km; stdout was:\n{out}')
    if quiet.strip():
        fails.append(f'(b) all stations matched but the reporter printed:\n{quiet}')

    if fails:
        print('FAIL test_offfault_station_dropped_report')
        for f in fails:
            print('  -', f)
        return 1
    print('SUCCESS test_offfault_station_dropped_report '
          f'(1 of 3 dropped -> named {named}; 0 dropped -> {len(quiet.strip())} bytes printed)')
    return 0


if __name__ == '__main__':
    sys.exit(main())

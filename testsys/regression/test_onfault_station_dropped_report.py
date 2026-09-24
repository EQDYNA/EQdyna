#! /usr/bin/env python3
"""
Regression guard (pathway_forward.md item 116) for ON-fault stations that
match no fault node being dropped with NO message -- the on-fault half of item
94 (test_offfault_station_dropped_report.py is the off-fault half and the
template this file follows).

THE BUG: `setOnFaultStation` (`src/fortran/meshgen.f90`) matches a station's
along-strike x and depth z against a fault node within tol, so a station off
the fault's node grid matches no node on any rank and `output_onfault_st`
writes no faultst* file for it -- silently. test.tpv29 at the 500 m gate dx
wrote 13 of 24 requested on-fault stations; test.meng2023a/cb at 400 m wrote
none of their 13 before PR #20 moved them onto the grid.

THE FIX: `checkOnFaultStationCoverage` (`eqdyna3d.f90`) OR-reduces the matched
station indices (`anonfs(2,:)`) over all ranks after meshgen, and rank 0 calls
`report_dropped_onfault_st` (`library_output.f90`), which prints a WARNING with
the dropped count and one line per dropped station (index, fault, x,z in km).


This test compiles `report_dropped_onfault_st` directly (no MPI, no case, no
mesh) and checks BEHAVIOUR, not source text (rule 10a):
  (a) 3 stations, station 2 unmatched -> the WARNING says "1 of 3", names
      station 2 (fault 1) with its x,z in km, and names neither 1 nor 3;
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
program item116_driver
    use globalvar
    implicit none
    logical :: matched(3)

    ntotft = 1
    allocate(nonfs(1), xonfs(2,3,1))
    nonfs(1) = 3
    xonfs(:,1,1) = (/     0.0d0,      0.0d0 /)
    xonfs(:,2,1) = (/ -18000.0d0, -15600.0d0 /)
    xonfs(:,3,1) = (/   5000.0d0, -12000.0d0 /)
    matched = (/ {m1}, {m2}, {m3} /)
    call report_dropped_onfault_st(matched, 3)
end program item116_driver
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
        print('FAIL test_onfault_station_dropped_report: gfortran not on PATH '
              '(this guard requires a real Fortran compile, it does not skip silently)')
        return 1
    tmp = tempfile.mkdtemp(prefix='testOnFaultDropped.')
    fails = []
    try:
        try:
            out = build_and_run(tmp, 'one', ('.true.', '.false.', '.true.'))
            quiet = build_and_run(tmp, 'none', ('.true.', '.true.', '.true.'))
        except RuntimeError as e:
            print('FAIL test_onfault_station_dropped_report')
            print(' -', e)
            return 1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    named = [int(n) for n in re.findall(r'dropped on-fault station (\d+)', out)]
    if not re.search(r'WARNING: 1 of 3 requested on-fault stations', out):
        fails.append(f'(a) no "WARNING: 1 of 3" line; stdout was:\n{out}')
    if named != [2]:
        fails.append(f'(a) named stations {named}, expected [2]; stdout was:\n{out}')
    m = re.search(r'station 2 \(fault 1\) at x,z =\s*(\S+)\s+(\S+)\s+km', out)
    if not m or [float(v) for v in m.groups()] != [-18.0, -15.6]:
        fails.append(f'(a) station 2 coordinates not reported as -18.000 -15.600 km; stdout was:\n{out}')
    if quiet.strip():
        fails.append(f'(b) all stations matched but the reporter printed:\n{quiet}')

    if fails:
        print('FAIL test_onfault_station_dropped_report')
        for f in fails:
            print('  -', f)
        return 1
    print('SUCCESS test_onfault_station_dropped_report '
          f'(1 of 3 dropped -> named {named}; 0 dropped -> {len(quiet.strip())} bytes printed)')
    return 0


if __name__ == '__main__':
    sys.exit(main())

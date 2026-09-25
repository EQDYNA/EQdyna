#! /usr/bin/env python3
"""
Regression guard (row 94 audit finding 4) for `output_offfault_st`
(`src/fortran/library_output.f90`): the header LOCATION STAMP must record
the ACTUAL matched node (`meshCoor(:,OffFaultStNodeIdIndex(2,i))`), not the
REQUESTED coordinate (`x4nds(:,OffFaultStNodeIdIndex(1,i))`) -- while the
FILE NAME stays keyed to the request (deliberately, see the fix's own
comment at library_output.f90:168-176).

THE GAP THIS CLOSES: `test_station_header_column_count.py`'s off-fault
driver sets `meshCoor(:,1) = x4nds(:,1)` (matched node == request), so a
revert of the header-stamp fix (back to stamping from x4nds) would still
pass that test -- the two sources are numerically identical there. This
test sets them to DIFFERENT values so the stamp's source is actually
distinguishable, and asserts on the ACTUAL, not the requested, value.

Cheap (rule 9): compiles 2 small .f90 files with plain gfortran once (no
MPI, no case, no mesh) and writes 1 tiny text file. Fails loudly (rule 2)
if gfortran is unavailable -- never silently skips.
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

# Requested station: x=500 m, y=-2000 m, z=-300 m (station 1's slot in
# OffFaultStNodeIdIndex is 1-indexed via column (1,i)).
# ACTUAL matched node (meshCoor, column OffFaultStNodeIdIndex(2,i)=7):
# x=600 m, y=-2100 m, z=-500 m -- deliberately DIFFERENT in every axis so
# the filename (request-derived) and the stamp (actual-derived) cannot be
# confused for one another by accident.
_DRIVER_SRC = r'''
program row94_finding4_offfault_driver
    use globalvar
    implicit none

    nstep = 1
    dx = 100.0d0
    dt = 0.01d0
    numOfOffFaultStCount = 1
    allocate(x4nds(3,1))
    x4nds(:,1) = (/ 500.0d0, -2000.0d0, -300.0d0 /)
    allocate(OffFaultStNodeIdIndex(2,1))
    OffFaultStNodeIdIndex(1,1) = 1
    OffFaultStNodeIdIndex(2,1) = 7
    allocate(meshCoor(3,7))
    meshCoor(:,7) = (/ 600.0d0, -2100.0d0, -500.0d0 /)
    allocate(OffFaultStGramSCEC(7,1))
    OffFaultStGramSCEC = 1.0d0

    call output_offfault_st
end program row94_finding4_offfault_driver
'''


def main():
    if shutil.which('gfortran') is None:
        print('FAIL test_offfault_station_header_actual_node: gfortran not on PATH '
              '(this guard requires a real Fortran compile, it does not skip silently)')
        return 1
    tmp = tempfile.mkdtemp(prefix='testOffFaultHeaderActualNode.')
    fails = []
    try:
        objs = []
        for fname in FORTRAN_DEPS:
            obj = os.path.join(tmp, fname.replace('.f90', '.o'))
            r = subprocess.run(
                ['gfortran', '-O0', '-ffree-line-length-none', '-I', tmp, '-J', tmp,
                 '-c', os.path.join(FSRC, fname), '-o', obj],
                capture_output=True, text=True, timeout=60)
            if r.returncode != 0:
                print('FAIL test_offfault_station_header_actual_node')
                print(' - compiling %s failed:\n%s\n%s' % (fname, r.stdout, r.stderr))
                return 1
            objs.append(obj)
        src = os.path.join(tmp, 'driver.f90')
        with open(src, 'w') as fh:
            fh.write(_DRIVER_SRC)
        binary = os.path.join(tmp, 'driver')
        r = subprocess.run(
            ['gfortran', '-O0', '-ffree-line-length-none', '-I', tmp, '-J', tmp]
            + objs + [src, '-o', binary], capture_output=True, text=True, timeout=60)
        if r.returncode != 0:
            print('FAIL test_offfault_station_header_actual_node')
            print(' - building driver failed:\n%s\n%s' % (r.stdout, r.stderr))
            return 1
        rundir = os.path.join(tmp, 'run')
        os.makedirs(rundir)
        r = subprocess.run([binary], cwd=rundir, capture_output=True, text=True, timeout=30)
        if r.returncode != 0:
            print('FAIL test_offfault_station_header_actual_node')
            print(' - driver exited %d:\n%s\n%s' % (r.returncode, r.stdout, r.stderr))
            return 1

        files = [f for f in os.listdir(rundir) if f.startswith('body')]
        if files != ['body-020st005dp003.txt']:
            fails.append('filename %r, expected [\'body-020st005dp003.txt\'] '
                         '(request-derived: y=-2000m->-20.0km->-020 (i4.3), x=500m->5.0km->005, '
                         'z=-300m->|3.0|km->003 -- NOT the actual node)' % files)
        else:
            with open(os.path.join(rundir, files[0])) as fh:
                text = fh.read()
            m = re.search(r'# location = (\S+) km off fault, (\S+) km along strike, (\S+) km depth', text)
            if not m:
                fails.append('no "# location = ... km off fault, ... km along strike, ... km '
                             'depth" line; file was:\n%s' % text)
            else:
                got = [float(v) for v in m.groups()]
                # Actual node: y=-2100m=-2.1km (off fault), x=600m=0.6km (along
                # strike), |z|=500m=0.5km (depth).
                want = [-2.1, 0.6, 0.5]
                if got != want:
                    fails.append('location stamp = %r, expected %r (the ACTUAL matched node, '
                                 'meshCoor column 7) -- got the REQUESTED (x4nds) coordinate '
                                 'instead if this shows -20.0/5.0/3.0' % (got, want))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    if fails:
        print('FAIL test_offfault_station_header_actual_node')
        for f in fails:
            print('  -', f)
        return 1
    print('SUCCESS test_offfault_station_header_actual_node '
          '(filename request-derived, header stamp actual-node-derived, and they differ)')
    return 0


if __name__ == '__main__':
    sys.exit(main())

#! /usr/bin/env python3
"""
Regression guard (pathway_forward.md item 85) for `output_onfault_st`
(`src/fortran/library_output.f90`) building a station-location header line
(`stLocStamp`, computed at :55) and never writing it.

THE BUG: `stLocStamp` was assigned every call and dropped -- no
`write(51,...)` of it existed anywhere in the subroutine. That silence had a
concrete cost: pathway item 67's own evidence command told the next reader to
find the on-fault station file "by its `# location = on fault` header line",
but no station file ever carried that line, so the command was unrunnable
from birth.

THE FIX emits the stamp as the header's first line, and -- because a value
that is built but never exercised has never been validated -- corrects the
down-dip field itself: `xonfs(2,...)` is vertical depth (matched against
`nodeCoor(3)` in meshgen.f90's `setOnFaultStation`, not a down-dip distance),
so the stamp must divide by sin(dip) to earn the word "down-dip", exactly as
the filename's own `dp` field two lines above it already does. Before the
fix this only went unnoticed because every case with on-fault stations that
is exercised elsewhere in this suite has dip=90 (sin=1) or C_degen==0
(fltxyz(2,4,*) forced to 90 regardless of true mesh dip); test.tpv36/
test.tpv37 (C_degen=15) have on-fault stations and a dip that is genuinely
not 90.

This test compiles `output_onfault_st` directly (same technique as
test_station_header_column_count.py) with a non-90 dip (30 degrees) and a
signed along-strike coordinate, then:

  1. Finds the line matching `# location = on fault, ... km along strike,
     ... km down-dip` -- the exact substring pathway item 67's evidence
     command greps for.
  2. Parses its two numbers and checks them against INDEPENDENTLY computed
     expected values (strike_m/1000, and depth_m/sin(dip)/1000 -- not
     depth_m/1000, which is the bug this guard would have caught) --
     pinning CONTENT, not merely presence, so a future refactor that emits
     `# location = ` with an empty or wrong value still fails.

Cheap (rule 9): compiles 2 small .f90 files with plain gfortran once (no
MPI, no case, no mesh) and writes 1 tiny text file; well under a second.
Fails loudly (rule 2) if gfortran is unavailable -- never silently skips.
"""
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FSRC = os.path.join(ROOT, 'src', 'fortran')

FORTRAN_DEPS = ('globalvar.f90', 'library_output.f90')

RE_LOCATION = re.compile(
    r'#\s*location\s*=\s*on fault,\s*([\-0-9.]+)\s*km along strike,\s*'
    r'([\-0-9.]+)\s*km down-dip')

# One station, one fault, one time step. dip=30 degrees (not 90) so the
# down-dip field must differ from the raw vertical-depth number to be
# correct; strike is negative to exercise the sign-aware field alongside it.
_DRIVER_SRC = r'''
program library_output_item85_driver
    use globalvar
    implicit none

    ntotft = 1
    friclaw = 1
    nstep = 1
    dx = 100.0d0
    dt = 0.01d0
    numOfOnFaultStCount = 1

    allocate(anonfs(3,1))
    anonfs(1,1) = 1; anonfs(2,1) = 1; anonfs(3,1) = 1

    allocate(xonfs(2,1,1))
    xonfs(1,1,1) = -2500.0d0; xonfs(2,1,1) = -3000.0d0

    allocate(fltxyz(2,4,1))
    fltxyz(2,4,1) = 30.0d0*pi/180.0d0

    allocate(onFaultQuantHistSCECForm(12,1,1))
    onFaultQuantHistSCECForm = 1.0d0

    call output_onfault_st
end program library_output_item85_driver
'''

STRIKE_M = -2500.0
DEPTH_M = -3000.0
DIP_DEG = 30.0

EXPECTED_STRIKE_KM = STRIKE_M / 1000.0
EXPECTED_DOWNDIP_KM = abs(DEPTH_M) / math.sin(math.radians(DIP_DEG)) / 1000.0
# What the pre-fix (uncorrected) computation would have produced -- used
# below to assert the guard actually distinguishes the two.
WRONG_DOWNDIP_KM = abs(DEPTH_M) / 1000.0


def build_and_run(tmp):
    """Compile globalvar.f90 + library_output.f90 + a tiny driver that calls
    output_onfault_st once, run it in `tmp`, and return the path of the one
    station file it wrote."""
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

    driver_src = os.path.join(tmp, 'driver.f90')
    with open(driver_src, 'w') as fh:
        fh.write(_DRIVER_SRC)
    driver_obj = os.path.join(tmp, 'driver.o')
    r = subprocess.run(
        ['gfortran', '-O0', '-ffree-line-length-none', '-I', tmp, '-J', tmp,
         '-c', driver_src, '-o', driver_obj],
        capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        raise RuntimeError(f'compiling driver failed:\n{r.stdout}\n{r.stderr}')

    binary = os.path.join(tmp, 'driver')
    r = subprocess.run(
        ['gfortran', '-O0'] + objs + [driver_obj, '-o', binary],
        capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        raise RuntimeError(f'linking driver failed:\n{r.stdout}\n{r.stderr}')

    rundir = os.path.join(tmp, 'run')
    os.makedirs(rundir, exist_ok=True)
    r = subprocess.run([binary], cwd=rundir, capture_output=True, text=True, timeout=30)
    if r.returncode != 0:
        raise RuntimeError(f'driver exited {r.returncode}:\n{r.stdout}\n{r.stderr}')

    files = [f for f in os.listdir(rundir) if f.startswith('faultst')]
    if len(files) != 1:
        raise RuntimeError(f'driver wrote {len(files)} faultst* files in {rundir}, expected 1')
    return os.path.join(rundir, files[0])


def main():
    if shutil.which('gfortran') is None:
        print('FAIL test_station_header_location_stamp: gfortran not on PATH '
              '(this guard requires a real Fortran compile, it does not skip silently)')
        return 1

    tmp = tempfile.mkdtemp(prefix='testStationHeaderLocationStamp.')
    try:
        try:
            path = build_and_run(tmp)
        except RuntimeError as e:
            print('FAIL test_station_header_location_stamp')
            print(' -', e)
            return 1

        with open(path) as fh:
            lines = fh.read().splitlines()

        match = None
        for line in lines:
            m = RE_LOCATION.search(line)
            if m:
                match = m
                break

        if match is None:
            print('FAIL test_station_header_location_stamp')
            print(f'  - no line matching "# location = on fault, ... km along '
                  f'strike, ... km down-dip" found in {path}')
            print('    (this is exactly the line pathway item 67\'s evidence '
                  'command greps for)')
            return 1

        got_strike_km = float(match.group(1))
        got_downdip_km = float(match.group(2))

        fails = []
        if abs(got_strike_km - EXPECTED_STRIKE_KM) > 0.05:
            fails.append(
                f'along-strike: expected {EXPECTED_STRIKE_KM:.1f} km, '
                f'found {got_strike_km:.1f} km (file: {path})')
        if abs(got_downdip_km - EXPECTED_DOWNDIP_KM) > 0.05:
            fails.append(
                f'down-dip: expected {EXPECTED_DOWNDIP_KM:.1f} km '
                f'(depth {abs(DEPTH_M)/1000.0:.1f} km / sin({DIP_DEG:.0f} deg)), '
                f'found {got_downdip_km:.1f} km (file: {path})')
        # Sanity: the two candidate values must actually differ, or this test
        # cannot tell a correct sin(dip) division from the pre-fix bug.
        if abs(EXPECTED_DOWNDIP_KM - WRONG_DOWNDIP_KM) < 0.5:
            fails.append(
                'test bug: EXPECTED_DOWNDIP_KM and WRONG_DOWNDIP_KM are too '
                'close to distinguish -- fixture is not exercising the '
                'sin(dip) correction')

        if fails:
            print('FAIL test_station_header_location_stamp')
            for f in fails:
                print('  -', f)
            return 1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print('SUCCESS test_station_header_location_stamp '
          f'(location line present; along-strike={EXPECTED_STRIKE_KM:.1f} km, '
          f'down-dip={EXPECTED_DOWNDIP_KM:.1f} km, both content-checked)')
    return 0


if __name__ == '__main__':
    sys.exit(main())

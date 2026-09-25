#! /usr/bin/env python3
"""
Regression guard (pathway_forward.md item 94) for off-fault station coverage:
a station whose requested depth is not a grid z-plane must still be WRITTEN,
at the nearest node, and named loudly -- not silently dropped.

THE ORIGINAL BUG: `setSurfaceStation` (`src/fortran/meshgen.f90`) matched a
station's depth EXACTLY (|nodeCoor(3) - x4nds(3,i)| < tol) and snapped only x
and y to the nearest interior node. A station whose depth was not a grid
z-plane (or that lay outside the interior x/y range) therefore matched no
node on any rank, and `output_offfault_st` wrote no file for it -- silently.
test.tpv8 at dx = 500 m requested 15 off-fault stations and wrote 11: its
four z = -0.3 km stations sit between the z = 0 and z = -0.5 km planes.

THE FIX (owner ruling 2026-09-24): depth now snaps to the nearest node the
same way x and y already did, so every station whose (x,y) is inside the
mesh matches SOME node. `report_dropped_offfault_st` (`library_output.f90`)
now distinguishes two cases:
  - matched, but not at the requested node: a SNAP -- "snapped off-fault
    station N (would otherwise be a dropped off-fault station): requested
    x,y,z = ... actual x,y,z = ... distance = ...".
  - matched on NO rank at all (x or y truly outside the mesh -- snapping
    depth cannot fix that): a true DROP -- "dropped off-fault station N at
    x,y,z = ...", unchanged wording from before this fix.

This test compiles `report_dropped_offfault_st` directly (no MPI, no case,
no mesh; same technique as test_station_header_column_count.py) and checks
BEHAVIOUR, not source text (rule 10a), with 3 stations:
  1: matched exactly at the requested node -> named nowhere.
  2: matched on no rank at all -> "WARNING: 1 of 3 ... match no grid node",
     named as "dropped off-fault station 2" with its requested km coordinates.
  3: matched, but 500 m off the requested depth (the clamp-to-nearest-node
     case this mission ports) -> "NOTICE: 1 of 3 ... do not sit exactly on a
     grid node", named as "snapped off-fault station 3 (would otherwise be a
     dropped off-fault station)" with requested vs actual (x,y,z) and a
     0.500 km distance.
  (b) all 3 matched at their exact requested node -> prints nothing at all.
Before this mission the subroutine took only `matchedAnyRank`, so the build
itself fails against the OLD signature. Fails loudly (rule 2) if gfortran is
unavailable -- never silently skips.
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

# Station 1: exact match (actual == requested) -> silent.
# Station 2: matched == .false. -> DROP. x4ndsZValidPersist(2) ({zvalid2})
#   selects which of the two DROP wordings (finding 6, row 94 audit):
#   .true. => depth was in-band, cause not checked here; .false. => the
#   CONFIRMED cause is the requested depth is outside the physical mesh band.
# Station 3: matched == .true. but actual z differs by 500 m -> SNAP.
_DRIVER_SRC = r'''
program item94_driver
    use globalvar
    implicit none
    logical :: matched(3)
    real (kind = dp) :: actual(3,3)

    totalNumOfOffSt = 3
    allocate(x4nds(3,3))
    allocate(x4ndsZValidPersist(3))
    x4nds(:,1) = (/    0.0d0,  1000.0d0,     0.0d0 /)
    x4nds(:,2) = (/    0.0d0,  -500.0d0,  -300.0d0 /)
    x4nds(:,3) = (/ 12000.0d0, 3000.0d0, -12000.0d0 /)
    x4ndsZValidPersist = (/ .true., {zvalid2}, .true. /)
    matched = (/ {m1}, {m2}, {m3} /)
    actual(:,1) = x4nds(:,1)
    ! actual(:,2) is never read when matched(2) is .false. (the dropped
    ! case) -- a real, sentinel-free value here just keeps the "all matched"
    ! run quiet without depending on that non-read to stay silent.
    actual(:,2) = x4nds(:,2)
    actual(:,3) = (/ 12000.0d0, 3000.0d0, {z3} /)
    call report_dropped_offfault_st(matched, actual)
end program item94_driver
'''


def build_and_run(tmp, tag, flags, z3, zvalid2='.true.'):
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
        fh.write(_DRIVER_SRC.format(m1=flags[0], m2=flags[1], m3=flags[2], z3=z3, zvalid2=zvalid2))
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
            out = build_and_run(tmp, 'one', ('.true.', '.false.', '.true.'), '-12500.0d0',
                                 zvalid2='.true.')
            quiet = build_and_run(tmp, 'none', ('.true.', '.true.', '.true.'), '-12000.0d0')
            # Finding 6 (row 94 audit): station 2 dropped with x4ndsZValidPersist
            # .false. -- the CONFIRMED cause (depth out of the physical band).
            confirmed = build_and_run(tmp, 'confirmed', ('.true.', '.false.', '.true.'),
                                       '-12500.0d0', zvalid2='.false.')
        except RuntimeError as e:
            print('FAIL test_offfault_station_dropped_report')
            print(' -', e)
            return 1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    dropped = [int(n) for n in re.findall(r'dropped off-fault station (\d+) at', out)]
    snapped = [int(n) for n in re.findall(r'snapped off-fault station (\d+)', out)]
    if not re.search(r'WARNING: 1 of 3 requested off-fault stations match no grid node', out):
        fails.append(f'no "WARNING: 1 of 3 ... match no grid node" line; stdout was:\n{out}')
    if dropped != [2]:
        fails.append(f'dropped stations {dropped}, expected [2]; stdout was:\n{out}')
    m = re.search(r'dropped off-fault station 2 at x,y,z =\s*(\S+)\s+(\S+)\s+(\S+)\s+km', out)
    if not m or [float(v) for v in m.groups()] != [0.0, -0.5, -0.3]:
        fails.append(f'station 2 coordinates not reported as 0.000 -0.500 -0.300 km; stdout was:\n{out}')
    if not re.search(r'NOTICE: 1 of 3 requested off-fault stations do not sit exactly on a grid z-plane', out):
        fails.append(f'no "NOTICE: 1 of 3 ... do not sit exactly on a grid z-plane" line; stdout was:\n{out}')
    if snapped != [3]:
        fails.append(f'snapped stations {snapped}, expected [3]; stdout was:\n{out}')
    if 'would otherwise be a dropped off-fault station' not in out:
        fails.append(f'snap line does not name itself a would-be drop; stdout was:\n{out}')
    m = re.search(r'snapped off-fault station 3 .*?requested x,y,z =\s*(\S+)\s+(\S+)\s+(\S+)\s+km, '
                  r'actual x,y,z =\s*(\S+)\s+(\S+)\s+(\S+)\s+km, distance =\s*(\S+)\s+km', out)
    if not m:
        fails.append(f'station 3 snap line not in the expected requested/actual/distance shape; stdout was:\n{out}')
    else:
        vals = [float(v) for v in m.groups()]
        if vals != [12.0, 3.0, -12.0, 12.0, 3.0, -12.5, 0.5]:
            fails.append(f'station 3 snap line values {vals}, expected '
                         f'[12.0, 3.0, -12.0, 12.0, 3.0, -12.5, 0.5]; stdout was:\n{out}')
    if 1 in dropped or 1 in snapped:
        fails.append(f'station 1 (exact match) named as dropped or snapped; stdout was:\n{out}')
    if quiet.strip():
        fails.append(f'all stations at their exact requested node but the reporter printed:\n{quiet}')

    # Finding 6 (row 94 audit): the DROP wording must distinguish a CHECKED
    # cause (depth outside the physical mesh band, x4ndsZValidPersist
    # .false.) from an unconfirmed one (x4ndsZValidPersist .true. -- depth
    # was fine, cause not further diagnosed), and must not claim "outside
    # the mesh" as if every drop's cause were checked.
    if 'match no grid node (outside the mesh)' in out or 'match no grid node (outside the mesh)' in confirmed:
        fails.append(f'summary line still claims every drop is confirmed outside the mesh; stdout:\n{out}\n{confirmed}')
    if 'cause not checked here' not in out:
        fails.append(f'station 2 (depth in-band, x4ndsZValidPersist=.true.) does not say its '
                     f'cause is unconfirmed; stdout was:\n{out}')
    if 'checked cause' in out:
        fails.append(f'station 2 (x4ndsZValidPersist=.true.) wrongly claims a checked cause; '
                     f'stdout was:\n{out}')
    if 'checked cause: requested depth is outside the physical, non-PML mesh band' not in confirmed:
        fails.append(f'station 2 (x4ndsZValidPersist=.false.) does not report the confirmed '
                     f'depth-out-of-band cause; stdout was:\n{confirmed}')
    if 'cause not checked here' in confirmed:
        fails.append(f'station 2 (x4ndsZValidPersist=.false.) wrongly says its cause is '
                     f'unconfirmed; stdout was:\n{confirmed}')

    if fails:
        print('FAIL test_offfault_station_dropped_report')
        for f in fails:
            print('  -', f)
        return 1
    print('SUCCESS test_offfault_station_dropped_report '
          f'(1 of 3 dropped -> named {dropped}; 1 of 3 snapped -> named {snapped}; '
          f'all exact -> {len(quiet.strip())} bytes printed)')
    return 0


if __name__ == '__main__':
    sys.exit(main())

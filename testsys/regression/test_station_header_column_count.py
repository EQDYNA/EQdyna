#! /usr/bin/env python3
"""
Regression guard (pathway_forward.md item 67) for `output_onfault_st`
(`src/fortran/library_output.f90`) declaring a column count that does not
match what it actually writes.

THE BUG: the '# Time series in 11 columns in format E15.7' line used to sit,
unconditional, ABOVE the `friclaw>=3` branch. The friclaw>=3 path does write
11 columns, but the friclaw<3 (slip-weakening) path at :95-108 writes only 8
-- t, h-slip, h-slip-rate, h-shear-stress, v-slip, v-slip-rate,
v-shear-stress, n-stress -- so every slip-weakening run declared 11 and
wrote 8. The 2015 SCEC archive (slip-weakening) declares 8 and writes 8.

THE FIX moved the declaration into each branch of the `if (friclaw>=3)`, so
it states that branch's own true count (11 or 8), and corrected the format
clause alongside it: the actual write is `E21.13` for column 1 and `E16.7`
for the rest in BOTH branches, never uniformly `E15.7` as the old string
claimed.

This test compiles `output_onfault_st` directly (same technique as
test_multifault_item9_fix.py: it depends on nothing but `globalvar`, no MPI,
no case, no mesh) once with `friclaw=1` and once with `friclaw=4`, then
checks THREE numbers against each other for the resulting station file,
exactly as a human would with:

    grep -o 'in [0-9]* columns' $f                          # declared count
    grep -A1 'names of the data fields' $f | tail -1 | wc -w # field-name count
    grep -v '^ *#' $f | head -1 | wc -w                      # data-value count

All three must agree (8/8/8 for friclaw<3, 11/11/11 for friclaw>=3). Before
the fix, the friclaw<3 run produced 11/8/8.

Cheap (rule 9): compiles 2 small .f90 files with plain gfortran twice (no
MPI, no case, no mesh) and writes 2 tiny text files; well under a second.
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

RE_DECLARED_NCOL = re.compile(r'[Tt]ime series in\s+(\d+)\s+columns')

# One station, one fault, one time step -- just enough for output_onfault_st
# to open unit 51, write its header, write one data row, and close it.
# `{friclaw}` and `{filename}` are substituted per driver build below;
# strike/dip differ between the two builds only so the two runs (in separate
# tmp dirs anyway) never risk producing the same filename by accident.
_DRIVER_SRC = r'''
program library_output_item67_driver
    use globalvar
    implicit none

    ntotft = 1
    friclaw = {friclaw}
    nstep = 1
    dx = 100.0d0
    dt = 0.01d0
    numOfOnFaultStCount = 1

    allocate(anonfs(3,1))
    anonfs(1,1) = 1; anonfs(2,1) = 1; anonfs(3,1) = 1

    allocate(xonfs(2,1,1))
    xonfs(1,1,1) = {strike_m}; xonfs(2,1,1) = -300.0d0

    allocate(fltxyz(2,4,1))
    fltxyz(2,4,1) = pi/2.0d0

    allocate(onFaultQuantHistSCECForm(12,1,1))
    onFaultQuantHistSCECForm = 1.0d0

    call output_onfault_st
end program library_output_item67_driver
'''


# Item 93 (2026-09-24): the OFF-fault writer, same three-number check, plus
# its filename. One station at y = 1.99 km: nint -> body020, the truncating
# int() the writer used to use -> body019.
_OFF_DRIVER_SRC = r'''
program library_output_item93_driver
    use globalvar
    implicit none

    nstep = 1
    dx = 100.0d0
    dt = 0.01d0
    numOfOffFaultStCount = 1
    allocate(x4nds(3,1))
    x4nds(:,1) = (/ 0.0d0, 1990.0d0, 0.0d0 /)
    allocate(OffFaultStNodeIdIndex(2,1))
    OffFaultStNodeIdIndex(1,1) = 1; OffFaultStNodeIdIndex(2,1) = 1
    ! Row 94: output_offfault_st's header stamp now reads meshCoor(:,
    ! OffFaultStNodeIdIndex(2,i)), the ACTUAL matched node -- exact match
    ! here (this test is about the column count/filename, not the stamp).
    allocate(meshCoor(3,1))
    meshCoor(:,1) = x4nds(:,1)
    allocate(OffFaultStGramSCEC(7,1))
    OffFaultStGramSCEC = 1.0d0

    call output_offfault_st
end program library_output_item93_driver
'''


def build_and_run_offfault(tmp):
    """Same build as build_and_run, driving output_offfault_st; returns the
    list of body* files it wrote."""
    objs = []
    for fname in FORTRAN_DEPS:
        obj = os.path.join(tmp, fname.replace('.f90', '.off.o'))
        r = subprocess.run(
            ['gfortran', '-O0', '-ffree-line-length-none', '-I', tmp, '-J', tmp,
             '-c', os.path.join(FSRC, fname), '-o', obj],
            capture_output=True, text=True, timeout=60)
        if r.returncode != 0:
            raise RuntimeError(f'compiling {fname} (off) failed:\n{r.stdout}\n{r.stderr}')
        objs.append(obj)
    src = os.path.join(tmp, 'driver.off.f90')
    with open(src, 'w') as fh:
        fh.write(_OFF_DRIVER_SRC)
    binary = os.path.join(tmp, 'driver.off')
    r = subprocess.run(['gfortran', '-O0', '-ffree-line-length-none', '-I', tmp, '-J', tmp]
                       + objs + [src, '-o', binary], capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        raise RuntimeError(f'building driver (off) failed:\n{r.stdout}\n{r.stderr}')
    rundir = os.path.join(tmp, 'run.off')
    os.makedirs(rundir, exist_ok=True)
    r = subprocess.run([binary], cwd=rundir, capture_output=True, text=True, timeout=30)
    if r.returncode != 0:
        raise RuntimeError(f'driver (off) exited {r.returncode}:\n{r.stdout}\n{r.stderr}')
    return rundir, sorted(f for f in os.listdir(rundir) if f.startswith('body'))


def check_offfault(tmp):
    fails = []
    rundir, files = build_and_run_offfault(tmp)
    if files != ['body020st000dp000.txt']:
        fails.append(f'off-fault: wrote {files}, expected [\'body020st000dp000.txt\'] '
                     f'(y=1.99 km must ROUND to 020 like the on-fault writer, not truncate to 019)')
    if len(files) != 1:
        return fails + [f'off-fault: expected exactly one body file, got {files}']
    got = three_counts(os.path.join(rundir, files[0]))
    if got != (7, 7, 7):
        fails.append(f'off-fault: declared/names/data = {got}, expected (7, 7, 7)')
    with open(os.path.join(rundir, files[0])) as fh:
        lines = fh.read().splitlines()
    names = [l for l in lines if l.strip().startswith('t h-disp')]
    if not names or not names[0].startswith(' t h-disp'):
        fails.append(f'off-fault: field-name line {names[:1]} is not the (1X,103A) shape '
                     f'the on-fault writer uses (one leading blank)')
    return fails


def build_and_run(tmp, friclaw, strike_m, tag):
    """Compile globalvar.f90 + library_output.f90 + a tiny driver that calls
    output_onfault_st once with the given friclaw, run it in `tmp`, and
    return the path of the one station file it wrote."""
    objs = []
    for fname in FORTRAN_DEPS:
        src = os.path.join(FSRC, fname)
        obj = os.path.join(tmp, fname.replace('.f90', f'.{tag}.o'))
        r = subprocess.run(
            ['gfortran', '-O0', '-ffree-line-length-none', '-I', tmp, '-J', tmp,
             '-c', src, '-o', obj],
            capture_output=True, text=True, timeout=60)
        if r.returncode != 0:
            raise RuntimeError(f'compiling {fname} ({tag}) failed:\n{r.stdout}\n{r.stderr}')
        objs.append(obj)

    driver_src = os.path.join(tmp, f'driver.{tag}.f90')
    with open(driver_src, 'w') as fh:
        fh.write(_DRIVER_SRC.format(friclaw=friclaw, strike_m=strike_m))
    driver_obj = os.path.join(tmp, f'driver.{tag}.o')
    r = subprocess.run(
        ['gfortran', '-O0', '-ffree-line-length-none', '-I', tmp, '-J', tmp,
         '-c', driver_src, '-o', driver_obj],
        capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        raise RuntimeError(f'compiling driver ({tag}) failed:\n{r.stdout}\n{r.stderr}')

    binary = os.path.join(tmp, f'driver.{tag}')
    r = subprocess.run(
        ['gfortran', '-O0'] + objs + [driver_obj, '-o', binary],
        capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        raise RuntimeError(f'linking driver ({tag}) failed:\n{r.stdout}\n{r.stderr}')

    rundir = os.path.join(tmp, f'run.{tag}')
    os.makedirs(rundir, exist_ok=True)
    r = subprocess.run([binary], cwd=rundir, capture_output=True, text=True, timeout=30)
    if r.returncode != 0:
        raise RuntimeError(f'driver ({tag}) exited {r.returncode}:\n{r.stdout}\n{r.stderr}')

    files = [f for f in os.listdir(rundir) if f.startswith('faultst')]
    if len(files) != 1:
        raise RuntimeError(f'driver ({tag}) wrote {len(files)} faultst* files in {rundir}, expected 1')
    return os.path.join(rundir, files[0])


def three_counts(path):
    """Reproduce, in Python, the three greps this guard is checking:
      grep -o 'in [0-9]* columns' $f
      grep -A1 'names of the data fields' $f | tail -1 | wc -w
      grep -v '^ *#' $f | head -1 | wc -w
    Returns (declared_ncol, name_count, data_count)."""
    with open(path) as fh:
        lines = fh.read().splitlines()

    declared = None
    for line in lines:
        m = RE_DECLARED_NCOL.search(line)
        if m:
            declared = int(m.group(1))
            break

    names_line = None
    for i, line in enumerate(lines):
        if 'names of the data fields' in line and i + 1 < len(lines):
            names_line = lines[i + 1]
            break
    name_count = len(names_line.split()) if names_line is not None else None

    data_line = None
    for line in lines:
        if not line.strip().startswith('#') and line.strip():
            data_line = line
            break
    data_count = len(data_line.split()) if data_line is not None else None

    return declared, name_count, data_count


def check_one(tmp, friclaw, strike_m, tag, expected):
    fails = []
    path = build_and_run(tmp, friclaw, strike_m, tag)
    declared, name_count, data_count = three_counts(path)
    got = (declared, name_count, data_count)
    if got != (expected, expected, expected):
        fails.append(
            f'friclaw={friclaw} ({tag}): declared={declared} names={name_count} '
            f'data={data_count}, expected all three == {expected} '
            f'(file: {path})')
    return fails


def main():
    if shutil.which('gfortran') is None:
        print('FAIL test_station_header_column_count: gfortran not on PATH '
              '(this guard requires a real Fortran compile, it does not skip silently)')
        return 1

    tmp = tempfile.mkdtemp(prefix='testStationHeaderColumnCount.')
    fails = []
    try:
        try:
            fails += check_one(tmp, friclaw=1, strike_m=500.0, tag='low', expected=8)
            fails += check_one(tmp, friclaw=4, strike_m=700.0, tag='high', expected=11)
            fails += check_offfault(tmp)
        except RuntimeError as e:
            print('FAIL test_station_header_column_count')
            print(' -', e)
            return 1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    if fails:
        print('FAIL test_station_header_column_count')
        for f in fails:
            print('  -', f)
        return 1
    print('SUCCESS test_station_header_column_count '
          '(friclaw<3 declares/names/data = 8/8/8, friclaw>=3 declares/names/data = 11/11/11, '
          'off-fault 7/7/7 in body020st000dp000.txt for y=1.99 km)')
    return 0


if __name__ == '__main__':
    sys.exit(main())

#! /usr/bin/env python3
"""
Regression guard (PROJECT_RULES.md rule 10; pathway_forward.md item 24(c)) for
the depth taper on the OFF-FAULT deviatoric pre-stress.

WHAT IT GUARDS. Until v5.9.0 `meshgen.f90`'s setPlasticStress built the
off-fault deviatoric pre-stress as a FIXED fraction of the vertical stress at
every depth (`devStr = abs(strVert)*devStrToStrVertRatio`), so the deviatoric
stress grew without bound with depth. SCEC TPV29/TPV30 (spec part 4, "Initial
Stress Tensor") instead taper the deviatoric component to zero between 17000 m
and 22000 m via

    Omega(depth) = 1                          , depth <= 17000 m
                 = (22000 - depth)/5000       , 17000 m <= depth <= 22000 m
                 = 0                          , depth >= 22000 m

`devStrDepthTaper` (src/fortran/func_lib.f90) is that Omega, with the two
depths as case input (par.devStrTaperDepthStart/End). This file checks it
three ways:

  1. Against the SPEC's closed form above -- the formula is transcribed from
     the spec document, not from the implementation.
  2. Against the Python port (`src/python/eqdyna/func_lib.dev_str_depth_taper`),
     for EXACT equality: both evaluate `(end - depth)/(end - start)` with the
     same operand order on the same doubles, so anything other than bitwise
     agreement is a real divergence, not rounding (CLAUDE.md: change physics in
     BOTH).
  3. That a case which does NOT configure a taper gets exactly 1.0 -- not
     0.999..., exactly 1.0. That is what makes the multiply in setPlasticStress
     an IEEE-exact identity and every already-gated case bit-for-bit unchanged.

HOW (real Fortran oracle, not a transcription): compiles a tiny driver against
the UNMODIFIED src/fortran/{globalvar,errorCodes,func_lib}.f90 and calls the
production `devStrDepthTaper` directly. Uses mpif90 because errorCodes.f90
(which func_lib.f90 uses) includes mpif.h; the driver itself never calls MPI,
so it runs as a plain serial process. Fails loudly (rule 2) if mpif90 is
absent -- it never silently skips.
"""
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FSRC = os.path.join(ROOT, 'src', 'fortran')
PYSRC = os.path.join(ROOT, 'src', 'python')
sys.path.insert(0, PYSRC)

FORTRAN_DEPS = ('globalvar.f90', 'errorCodes.f90', 'func_lib.f90')

_DRIVER_SRC = r'''
program devStrTaperDriver
    use globalvar
    implicit none
    real (kind = dp) :: devStrDepthTaper
    real (kind = dp) :: depth
    integer :: n, i

    read(*,*) devStrTaperDepthStart, devStrTaperDepthEnd
    read(*,*) n
    do i = 1, n
        read(*,*) depth
        write(*,'(ES25.16E3)') devStrDepthTaper(depth)
    enddo
end program devStrTaperDriver
'''

# TPV29/TPV30's own taper depths (spec part 4), and the "not configured" pair
# scripts/case.setup writes for every other case.
TPV30_TAPER = (17000.0, 22000.0)
NO_TAPER = (0.0, 0.0)

# Depths (m, positive down) probed in both configurations: above the taper,
# either side of and exactly on both knees, inside the ramp, below the taper,
# and one above the free surface (a negative depth must not produce a
# >1 multiplier).
DEPTHS = (-100.0, 0.0, 8000.0, 16999.9, 17000.0, 17000.1, 19500.0,
          21999.9, 22000.0, 22000.1, 25000.0, 40000.0)


def specOmega(depth, start, end):
    """SCEC TPV29/30's Omega(depth), transcribed from the spec text, for the
    general (start, end) the code takes; `end <= start` means no taper."""
    if not end > start:
        return 1.0
    if depth <= start:
        return 1.0
    if depth >= end:
        return 0.0
    return (end - depth)/(end - start)


def buildDriver(tmp):
    objs = []
    for fname in FORTRAN_DEPS:
        obj = os.path.join(tmp, fname.replace('.f90', '.o'))
        r = subprocess.run(
            ['mpif90', '-O0', '-ffree-line-length-none', '-I', tmp, '-J', tmp,
             '-c', os.path.join(FSRC, fname), '-o', obj],
            capture_output=True, text=True, timeout=120)
        if r.returncode != 0:
            raise RuntimeError(f'compiling {fname} failed:\n{r.stdout}\n{r.stderr}')
        objs.append(obj)

    src = os.path.join(tmp, 'devStrTaperDriver.f90')
    with open(src, 'w') as fh:
        fh.write(_DRIVER_SRC)
    obj = os.path.join(tmp, 'devStrTaperDriver.o')
    r = subprocess.run(
        ['mpif90', '-O0', '-ffree-line-length-none', '-I', tmp, '-J', tmp,
         '-c', src, '-o', obj],
        capture_output=True, text=True, timeout=120)
    if r.returncode != 0:
        raise RuntimeError(f'compiling the driver failed:\n{r.stdout}\n{r.stderr}')

    binary = os.path.join(tmp, 'devStrTaperDriver')
    r = subprocess.run(['mpif90', '-O0'] + objs + [obj, '-o', binary],
                       capture_output=True, text=True, timeout=120)
    if r.returncode != 0:
        raise RuntimeError(f'linking the driver failed:\n{r.stdout}\n{r.stderr}')
    return binary


def runFortran(binary, taper, depths):
    stdin = '%r %r\n%d\n' % (taper[0], taper[1], len(depths))
    stdin += ''.join('%r\n' % d for d in depths)
    r = subprocess.run([binary], input=stdin, capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        raise RuntimeError(f'driver exited {r.returncode}:\n{r.stdout}\n{r.stderr}')
    vals = [float(ln) for ln in r.stdout.split()]
    if len(vals) != len(depths):
        raise RuntimeError(f'expected {len(depths)} values, got {vals!r}')
    return vals


def main():
    if shutil.which('mpif90') is None:
        print('FAIL test_dev_str_depth_taper: mpif90 not on PATH '
              '(this guard requires a real Fortran compile, it does not skip silently)')
        return 1

    import numpy as np
    from eqdyna.func_lib import dev_str_depth_taper

    print('Regression guard: off-fault deviatoric pre-stress depth taper')
    fails = []
    tmp = tempfile.mkdtemp(prefix='testDevStrTaper.')
    try:
        try:
            binary = buildDriver(tmp)
        except RuntimeError as e:
            print('FAIL test_dev_str_depth_taper')
            print(' -', e)
            return 1

        for label, taper in (('no taper (case.setup default)', NO_TAPER),
                              ('TPV30 17-22 km', TPV30_TAPER)):
            fort = runFortran(binary, taper, DEPTHS)
            port = [float(v) for v in
                    dev_str_depth_taper(np.asarray(DEPTHS), taper[0], taper[1])]
            worstSpec = 0.0
            for depth, f, p in zip(DEPTHS, fort, port):
                expected = specOmega(depth, taper[0], taper[1])
                worstSpec = max(worstSpec, abs(f - expected))
                if abs(f - expected) > 1e-15:
                    fails.append(f'{label}: depth {depth} m -- fortran {f!r} vs '
                                  f'spec Omega {expected!r}')
                if f != p:
                    fails.append(f'{label}: depth {depth} m -- fortran {f!r} vs '
                                  f'python port {p!r} (must be bitwise equal)')
            print(f'  {label}: {len(DEPTHS)} depths, max |fortran - spec| = '
                  f'{worstSpec:.3e}, fortran == port bitwise')

        # The identity that makes every already-gated case bit-for-bit: an
        # unconfigured taper must return exactly 1.0, at every depth.
        offFort = runFortran(binary, NO_TAPER, DEPTHS)
        if any(v != 1.0 for v in offFort):
            fails.append(f'no taper: expected exactly 1.0 at every depth, got {offFort!r}')
        else:
            print('  no taper returns exactly 1.0 (IEEE-exact multiplicative identity)')

        # A configured taper must actually change something inside the ramp --
        # otherwise the guard above would pass on a no-op implementation.
        onFort = runFortran(binary, TPV30_TAPER, (19500.0,))
        if onFort[0] != 0.5:
            fails.append('TPV30 taper: Omega at the ramp midpoint 19500 m should be '
                          f'exactly 0.5, got {onFort[0]!r}')
        else:
            print('  TPV30 taper bites: Omega(19500 m) == 0.5 exactly')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    if fails:
        print('FAIL test_dev_str_depth_taper')
        for f in fails:
            print('  -', f)
        return 1
    print('SUCCESS test_dev_str_depth_taper')
    return 0


if __name__ == '__main__':
    sys.exit(main())

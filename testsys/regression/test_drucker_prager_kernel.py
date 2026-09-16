#! /usr/bin/env python3
"""
Regression guard (PROJECT_RULES.md rule 10; pathway_forward.md item 23) for
the Drucker-Prager viscoplastic return-mapping kernel having ZERO isolated
test coverage.

Background: `calcElemKU.f90:127-161` (the `if (C_elastic==0)` block) is
exercised today only end-to-end, inside test.drv.a6 and test.tpv104, which
sit behind a chaos-tolerant flip-budget gate (see
testsys/parity/evidence_drv_a6_chaos.py / matrix.DRV_A6). That gate tolerates
up to 450 rupture-arrival flips out of 5151 fault nodes as ordinary
discretisation/summation-order noise -- so a genuine math error in the D-P
block (wrong sign, swapped operand, mis-clamped yield, ...) could easily land
INSIDE that budget and ship silently. This item was flagged during the
2026-09-14 release audit (item 21) and never closed with an actual test
(item 23) -- this file closes it.

WHAT THIS TESTS AND HOW (real Fortran oracle, not a transcription):
This compiles a tiny standalone driver (`_DRIVER_SRC` below) against the
UNMODIFIED `src/fortran/{globalvar,calcB,calcQAttenuationCoeff,calcElemKU}.f90`
and calls the REAL, production `calcElemKU` subroutine directly -- it does
not re-type or reimplement any of the D-P arithmetic in a new Fortran
snippet. The driver's only job is to give `calcElemKU` inputs that make every
other physics path in that subroutine a no-op:
  - vl = dl = 0 (nodal velocity/displacement) => strainrate = strain = 0
    => `stress(i) = stress(i) + stressrate(i)*dt` (calcElemKU.f90:74) leaves
    stress(1:6) EXACTLY as the caller set it, for any dt/mate.
  - C_Q = 0 selects the plain (non-attenuating) stress-update branch.
  - C_elastic = 0 is what turns the D-P block on at all.
So the ONLY thing that can move stress(1:6)/pstrmag away from the values
fed in on stdin is `calcElemKU.f90:127-161` itself -- the D-P block is
isolated by construction, not by extraction. Verified: with these inputs and
a stress state built to sit below yield, the driver returns stress
UNCHANGED and pstrmag == 0 (see CASES['elastic'] below) -- i.e. calling the
real subroutine this way really does start from "no-op" and let only the
plastic branch move anything.

This is compared against the Python port's `_drucker_prager`
(src/python/eqdyna/assembleGlobalKU.py:307), called directly and in
isolation (no solver loop, no mesh, no case) with the SAME synthetic inputs.

OUT OF SCOPE, DOCUMENTED, NOT A BUG: the Python port never computes
pstrain/pstrmag (assembleGlobalKU.py:316-319 explains why: it is write-only
in the Fortran and unverifiable against any committed reference). This test
prints the Fortran pstrmag for diagnostic purposes and uses it ONLY to check
that each case actually exercises the branch its name claims (elastic cases
must show pstrmag==0, yielding cases must show pstrmag>0) -- it is never
compared against the Python port.

TOLERANCE: rtol=1e-12, atol=1e-9 (Pa; ~1e-16 relative to the 1e6-1e7 Pa
stresses used here). This is deterministic scalar IEEE-754 arithmetic
end-to-end -- taomax's second-invariant sum and yld's clamp use the IDENTICAL
operand order in both languages, and IEEE-754 mandates correctly-rounded
sqrt in any conformant implementation, so those contribute exactly zero
divergence. exp() is the one place bit-identity is not guaranteed by the
standard (gfortran's libm call vs NumPy's own vectorized exp), which is why
this uses a tight-but-not-bitwise tolerance rather than `==`.

Cheap (rule 9): compiles 4 small non-MPI .f90 files with plain gfortran
(no mpif90, no case, no mesh) and runs a handful of scalar cases -- well
under a second total. Fails loudly (rule 2) if gfortran is unavailable --
never silently skips.
"""
import os
import subprocess
import shutil
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FSRC = os.path.join(ROOT, 'src', 'fortran')
PYSRC = os.path.join(ROOT, 'src', 'python')

sys.path.insert(0, PYSRC)

# Fortran source files the driver links against, UNMODIFIED, straight out of
# the real solver tree -- calcElemKU.f90 is the file under test, the other
# three are what it needs to compile/link (calcB for the strain-displacement
# matrix it never uses meaningfully here since vl=dl=0, calcQAttenuationCoeff
# only because calcElemKU references it in a branch this driver never takes
# but the linker still needs the symbol).
FORTRAN_DEPS = ('globalvar.f90', 'calcB.f90', 'calcQAttenuationCoeff.f90', 'calcElemKU.f90')

# The driver: sets the handful of globalvar module scalars calcElemKU.f90's
# D-P block reads (dt, ccosphi, sinphi, tv), feeds it a caller-chosen
# stress(1:6) with vl=dl=0/C_Q=0 so nothing else in the subroutine can touch
# it, and prints stress(1:6) + pstrmag at full double precision. One process,
# one stdin/stdout protocol, no case, no mesh, no MPI.
_DRIVER_SRC = r'''
program dp_kernel_driver
    implicit none
    call run_cases()
end program dp_kernel_driver

subroutine run_cases()
    use globalvar
    implicit none
    integer :: ncases, ic, i
    real(kind=dp) :: globalShapeFunc(nrowsh-1, nen)
    real(kind=dp) :: mate(5), vl(nee), dl(nee), stress(12), elresf(nee)
    real(kind=dp) :: constk, porep, pstrmag, ex(3,8)
    real(kind=dp) :: s_dt, s_ccosphi, s_sinphi, s_tv
    real(kind=dp) :: s1, s2, s3, s4, s5, s6

    C_Q = 0
    C_elastic = 0
    w = 0.0d0
    rdampk = 0.0d0
    globalShapeFunc = 0.0d0
    vl = 0.0d0
    dl = 0.0d0
    ex = 0.0d0
    mate = 0.0d0
    mate(4) = 3.0d10   ! lam: multiplies strainrate==0 only, never reaches output
    mate(5) = 3.0d10   ! miu: only scales pstrmag, which this test never gates on
    constk = 1.0d0
    porep = 0.0d0       ! matches the port's hardcoded 0.0 (assembleGlobalKU.py:311-314)

    read(*,*) ncases
    do ic = 1, ncases
        read(*,*) s_dt, s_ccosphi, s_sinphi, s_tv, s1, s2, s3, s4, s5, s6
        dt = s_dt
        ccosphi = s_ccosphi
        sinphi = s_sinphi
        tv = s_tv
        stress = 0.0d0
        stress(1) = s1; stress(2) = s2; stress(3) = s3
        stress(4) = s4; stress(5) = s5; stress(6) = s6
        elresf = 0.0d0
        pstrmag = 0.0d0

        call calcElemKU(globalShapeFunc, mate, vl, dl, stress, elresf, constk, porep, pstrmag, ex)

        write(*,'(7ES25.16E3)') (stress(i), i=1,6), pstrmag
    end do
end subroutine run_cases
'''

# Synthetic (dt, ccosphi, sinphi, tv, stress(1:6)) states, one per named
# regime. `expect_yield` is a self-check on the case DESIGN (asserted against
# the real Fortran pstrmag below), not a hand-computed expected numeric
# answer -- the numeric answer comes only from the compiled oracle.
CASES = {
    'elastic': dict(
        dt=0.01, ccosphi=5.0e6, sinphi=0.6, tv=0.05,
        stress=(-1.0e7, -1.0e7, -1.0e7, 0.0, 0.0, 0.0),
        expect_yield=False,
        note='isotropic compression, taomax=0 << yield=1.1e7: no adjustment',
    ),
    'yielding': dict(
        dt=0.01, ccosphi=2.0e6, sinphi=0.6, tv=0.05,
        stress=(-6.0e7, -2.0e7, 2.0e7, 0.0, 0.0, 0.0),
        expect_yield=True,
        note='taomax=4.0e7 >> yield=1.4e7: well inside the plastic branch',
    ),
    'near_boundary': dict(
        dt=0.02, ccosphi=3.0e6, sinphi=0.5, tv=0.1,
        stress=(-4499989.5, -15000000.0, -25500010.5, 0.0, 0.0, 0.0),
        expect_yield=True,
        note='taomax=yield*(1+1e-6): exercises the taomax>yield branch edge, '
             'not the exactly-equal case (which is not well-defined to test)',
    ),
    'shear_yielding': dict(
        dt=0.005, ccosphi=1.0e6, sinphi=0.4, tv=0.02,
        stress=(-1.5e7, -1.5e7, -1.5e7, 8.0e6, -5.0e6, 3.0e6),
        expect_yield=True,
        note='strdev(1:3)=0, taomax carried entirely by strdev(4:6): exercises '
             'the shear terms of the second-invariant formula, not just normal stress',
    ),
}
CASE_ORDER = ('elastic', 'yielding', 'near_boundary', 'shear_yielding')

RTOL = 1e-12
ATOL = 1e-9


def build_driver(tmp):
    """Compile the driver against the REAL, unmodified calcElemKU.f90 (and
    the three files it needs to link) with plain gfortran -- no MPI, no
    case, no mesh. Raises loudly on any compiler/linker failure."""
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

    driver_src = os.path.join(tmp, 'dp_kernel_driver.f90')
    with open(driver_src, 'w') as fh:
        fh.write(_DRIVER_SRC)
    driver_obj = os.path.join(tmp, 'dp_kernel_driver.o')
    r = subprocess.run(
        ['gfortran', '-O0', '-ffree-line-length-none', '-I', tmp, '-J', tmp,
         '-c', driver_src, '-o', driver_obj],
        capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        raise RuntimeError(f'compiling driver failed:\n{r.stdout}\n{r.stderr}')

    binary = os.path.join(tmp, 'dp_kernel_driver')
    r = subprocess.run(
        ['gfortran', '-O0'] + objs + [driver_obj, '-o', binary],
        capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        raise RuntimeError(f'linking driver failed:\n{r.stdout}\n{r.stderr}')
    return binary


def run_fortran(binary):
    """Feed all CASES to the compiled driver in one process; return
    {case_name: (stress1..stress6, pstrmag)}."""
    lines = [str(len(CASE_ORDER))]
    for name in CASE_ORDER:
        c = CASES[name]
        s = c['stress']
        lines.append(f"{c['dt']!r} {c['ccosphi']!r} {c['sinphi']!r} {c['tv']!r} "
                      f"{s[0]!r} {s[1]!r} {s[2]!r} {s[3]!r} {s[4]!r} {s[5]!r}")
    stdin = '\n'.join(lines) + '\n'
    r = subprocess.run([binary], input=stdin, capture_output=True, text=True, timeout=30)
    if r.returncode != 0:
        raise RuntimeError(f'driver exited {r.returncode}:\n{r.stdout}\n{r.stderr}')
    out_lines = [ln for ln in r.stdout.strip().splitlines() if ln.strip()]
    if len(out_lines) != len(CASE_ORDER):
        raise RuntimeError(f'expected {len(CASE_ORDER)} output lines, got '
                            f'{len(out_lines)}:\n{r.stdout}')
    results = {}
    for name, line in zip(CASE_ORDER, out_lines):
        vals = [float(x) for x in line.split()]
        if len(vals) != 7:
            raise RuntimeError(f'case {name}: expected 7 values, got {vals!r}')
        results[name] = (tuple(vals[:6]), vals[6])
    return results


def run_python(name):
    """Call the port's _drucker_prager directly (no solver loop) on the same
    synthetic inputs, and return the resulting stress(1:6) as a 6-tuple."""
    import numpy as np
    from eqdyna.assembleGlobalKU import _drucker_prager

    c = CASES[name]
    inv = {'ccosphi': c['ccosphi'], 'sinphi': c['sinphi'], 'tv': c['tv']}
    stress_i = np.array([c['stress']], dtype=np.float64)
    out = _drucker_prager(np, inv, stress_i, c['dt'])
    return tuple(float(v) for v in out[0])


def main():
    if shutil.which('gfortran') is None:
        print('FAIL test_drucker_prager_kernel: gfortran not on PATH '
              '(this guard requires a real Fortran compile, it does not skip silently)')
        return 1

    fails = []
    tmp = tempfile.mkdtemp(prefix='testDrPragerKernel.')
    try:
        try:
            binary = build_driver(tmp)
        except RuntimeError as e:
            print('FAIL test_drucker_prager_kernel')
            print(' -', e)
            return 1

        try:
            fortran_results = run_fortran(binary)
        except RuntimeError as e:
            print('FAIL test_drucker_prager_kernel')
            print(' -', e)
            return 1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    for name in CASE_ORDER:
        c = CASES[name]
        fort_stress, fort_pstrmag = fortran_results[name]

        # Self-check on the CASE DESIGN (not on the port): confirm each named
        # regime actually exercised the branch its name claims, per the real
        # Fortran oracle's own pstrmag.
        yielded = fort_pstrmag > 0.0
        if yielded != c['expect_yield']:
            fails.append(
                f'{name}: case design check failed -- expected yield={c["expect_yield"]} '
                f'but Fortran oracle reports pstrmag={fort_pstrmag!r} ({c["note"]})')
            continue

        port_stress = run_python(name)
        for k in range(6):
            f, p = fort_stress[k], port_stress[k]
            tol = ATOL + RTOL * abs(f)
            if abs(f - p) > tol:
                fails.append(
                    f'{name}: stress({k + 1}) mismatch -- fortran={f!r} python={p!r} '
                    f'abs_diff={abs(f - p)!r} (tol={tol!r}); {c["note"]}')

    if fails:
        print('FAIL test_drucker_prager_kernel')
        for f in fails:
            print('  -', f)
        return 1
    print(f'SUCCESS test_drucker_prager_kernel ({len(CASE_ORDER)} cases, '
          f'fortran oracle == python port to rtol={RTOL:g}/atol={ATOL:g})')
    return 0


if __name__ == '__main__':
    sys.exit(main())


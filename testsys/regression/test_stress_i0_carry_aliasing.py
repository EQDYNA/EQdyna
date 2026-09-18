#! /usr/bin/env python3
"""
Regression guard (pathway_forward.md item 41) for a latent aliasing bug in
driver.py's time-loop carry construction.

`driver.py`'s `run()` seeds the per-step carry state with:

    xp.asarray(inv['stress_i0'])

Under the numpy backend, `np.asarray` on an argument that is ALREADY an
ndarray returns the SAME object -- no copy. `inv['stress_i0']` is built once
in `assembleGlobalKU.py` (`build()`, around init_stress[E_int].copy()) and,
before this fix, that one line was the array's only read site anywhere in the
shipped codebase, so the aliasing was harmless in practice. It is latent, not
inert: `driver.py` mutates the carry slot IN PLACE at
`assembleGlobalKU.py:261` (deliberate, for performance) every step, so a
second reader of `inv['stress_i0']` -- present or future -- would silently
observe corrupted, post-simulation state instead of the mesh's t=0 initial
stress.

This test proves the fix (`xp.asarray(inv['stress_i0']).copy()`) two ways,
against the REAL production code path (`driver.run`), not a re-typed stand-in:
  1. object identity: the array `driver.run` assigns into the carry tuple is
     NOT `inv['stress_i0']`.
  2. isolation: mutating the carry array in place (numpy backend -- jax
     arrays are immutable, so this half is a no-op there by construction) does
     NOT change `inv['stress_i0']`.

HOW, without paying for a real time-stepping run: `backend.run_time_loop` is
monkeypatched to capture `(inv, carry0)` and raise a sentinel exception BEFORE
any stepping happens -- `driver.run` builds the real invariants and the real
carry tuple exactly as production does, and this only skips the (expensive,
irrelevant-to-this-bug) actual physics loop.

Cheap (rule 9): one tiny serial test.tpv8 case, pure Python (no Fortran build,
no MPI) -- create.newcase + case.setup + build_solver_state, ~2 s measured on
this box; no actual time-stepping happens for either backend.
"""
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PYSRC = os.path.join(ROOT, 'src', 'python')
sys.path.insert(0, PYSRC)

CASE = 'test.tpv8'


class _StopBeforeStepping(Exception):
    """Raised by the monkeypatched run_time_loop to bail out after the real
    invariants/carry state are built, before any (expensive, irrelevant-to-
    this-bug) actual time-stepping runs."""


def _env():
    env = dict(os.environ)
    env['EQDYNAROOT'] = ROOT
    env['PATH'] = os.pathsep.join([os.path.join(ROOT, 'bin'),
                                   os.path.join(ROOT, 'scripts'),
                                   env.get('PATH', '')])
    return env


def _make_serial_case(case_dir):
    """create.newcase + force a SERIAL decomposition + case.setup, same
    pattern as testsys/e2e/run_e2e.py:make_serial_case -- no Fortran binary
    involved, case.setup writes the bFile/netCDF inputs the standalone
    Python solver reads natively."""
    env = _env()
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'create.newcase'),
                        case_dir, CASE], env=env, capture_output=True, text=True)
    if r.returncode != 0:
        raise AssertionError('create.newcase failed: %s' % (r.stderr or '')[-800:])
    params = os.path.join(case_dir, 'user_defined_params.py')
    with open(params) as f:
        text = f.read().rstrip('\n')
    with open(params, 'w') as f:
        f.write(text + '\n\n# forced serial by test_stress_i0_carry_aliasing.py'
                       ' (the standalone Python solver is serial-only)\n'
                       'par.nx = 1\npar.ny = 1\npar.nz = 1\n')
    r = subprocess.run([sys.executable, 'case.setup'], cwd=case_dir,
                       env=env, capture_output=True, text=True)
    if r.returncode != 0:
        raise AssertionError('case.setup failed: %s' % (r.stderr or '')[-800:])


def _capture_inv_and_carry0(S, xp):
    """Run the REAL driver.run() far enough to build inv/carry0 exactly as
    production does, then intercept them via a monkeypatched
    backend.run_time_loop that raises before any stepping. Returns
    (inv, carry0)."""
    from eqdyna import backend as B
    from eqdyna import driver

    captured = {}
    real_run_time_loop = B.run_time_loop

    def fake_run_time_loop(xp_, build_step, inv, carry, n):
        captured['inv'] = inv
        captured['carry0'] = carry
        raise _StopBeforeStepping()

    B.run_time_loop = fake_run_time_loop
    try:
        try:
            driver.run(S, nsteps=1, verbose=False, xp=xp)
        except _StopBeforeStepping:
            pass
        else:
            raise AssertionError('driver.run() finished without reaching '
                                  'run_time_loop -- the interception point moved')
    finally:
        B.run_time_loop = real_run_time_loop

    if 'inv' not in captured or 'carry0' not in captured:
        raise AssertionError('run_time_loop was never called by driver.run()')
    return captured['inv'], captured['carry0']


# carry0's element-4 slot is `xp.asarray(inv['stress_i0']).copy()` --
# driver.py's carry0 tuple order: (v1, velArr, dispArr, force, stress_i,
# s_p, fric, fnft, timeElapsed, ...).
STRESS_I_SLOT = 4


def _check_backend(backend_name):
    from eqdyna import backend as B
    xp = B.array_module(backend_name)

    with tempfile.TemporaryDirectory() as tmp:
        case_dir = os.path.join(tmp, 'case')
        _make_serial_case(case_dir)

        from eqdyna import eqdyna3d as E
        S, _mesh = E.build_solver_state(case_dir)

        inv, carry0 = _capture_inv_and_carry0(S, xp)
        carry_stress = carry0[STRESS_I_SLOT]
        inv_stress = inv['stress_i0']

        if carry_stress is inv_stress:
            raise AssertionError(
                '%s: carry0[%d] IS inv["stress_i0"] (same object) -- the '
                'asarray-without-copy aliasing bug is back' % (backend_name, STRESS_I_SLOT))

        if backend_name == 'numpy':
            import numpy as np
            before = inv_stress.copy()
            carry_stress += 1.0   # the exact style of in-place mutation
                                  # assembleGlobalKU.py:261 performs every step
            if not np.array_equal(inv_stress, before):
                raise AssertionError(
                    'numpy: mutating the carry stress array in place changed '
                    'inv["stress_i0"] -- aliasing bug reproduced')
            if np.array_equal(carry_stress, before):
                raise AssertionError(
                    'numpy: mutation of the carry array had no effect -- test '
                    'is not actually exercising in-place mutation')

    print('  %s: carry0[%d] is a distinct object from inv["stress_i0"]%s'
          % (backend_name, STRESS_I_SLOT,
             ', and in-place mutation does not leak into it' if backend_name == 'numpy' else ''))


def main():
    print('Regression guard: driver.py stress_i0 carry-state aliasing (pathway item 41)')
    fails = []
    for backend_name in ('numpy', 'jax'):
        try:
            _check_backend(backend_name)
        except Exception as e:
            fails.append('%s: %s' % (backend_name, e))

    if fails:
        print('FAIL test_stress_i0_carry_aliasing')
        for f in fails:
            print('  -', f)
        return 1
    print('SUCCESS test_stress_i0_carry_aliasing (numpy and jax backends)')
    return 0


if __name__ == '__main__':
    sys.exit(main())

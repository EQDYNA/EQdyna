#! /usr/bin/env python3
"""
Standalone-Python-EQdyna END-TO-END acceptance tier (`testsys/run.py accept`).

Makes the tpv8/tpv104/tpv1053d acceptance numbers previously recorded ONLY
in a commit message (M7.5/M8 landing commit 057e589; PROJECT_RULES rule 2
"no silent placeholder"/rule 4 "only fresh runs are evidence" -- a number
that only a human can re-derive by reading git log is neither) into a
reproducible, gated, re-runnable check.

For each of test.tpv8 (friclaw=1), test.tpv104 (friclaw=4), test.tpv1053d
(friclaw=5), IN SEQUENCE (never concurrently -- this is a shared machine,
same constraint as testsys/parity/run_parity.py and make_fixtures.py):

  (a) `create.newcase` + `case.setup` a SERIAL (nx=ny=nz=1) case in a
      fresh temp directory (no Fortran binary needed here -- the
      standalone solver reads only the case-input FILES case.setup
      writes, never runs eqdyna).
  (b) `python3 -m eqdyna.standalone <case_dir>` with the DEFAULT (jax)
      backend -- the same invocation a real user would type.
  (c) Coordinate-aligned comparison against the COMMITTED
      test.reference.results/<case>/frt.txt0 + frt.txt2 (a 4-rank
      reference; a serial run's frt.txt0 uses a different domain
      decomposition and therefore different node order/count, so a flat
      positional diff cannot compare the two directly -- see
      python/eqdyna/standalone/main.py's module docstring): dedupe the
      combined reference by rounded (x,y,z) coordinates (partition-
      boundary nodes are written by both ranks), lexsort both the
      deduped reference and the python output by (x,y,z), then
      np.allclose at PROJECT_RULES's one calibrated threshold (1e-3,
      same as check.test.py -- rule 5, never an invented tolerance).
  (d) prints `SUCCESS <case> max_abs_diff=<n>` or `FAIL <case> ...` and
      accumulates a non-zero exit if ANY case fails -- no case is allowed
      to fail silently or be skipped without being reported as a failure.

Verified fresh in the session that authored this file: tpv8 8.23e-11,
tpv104 3.39e-7, tpv1053d 6.65e-6 (the last of these AFTER discovering and
fixing a stale `bin/eqdyna` that had pre-dated the fric_tp_h fix -- see
pathway_forward.md/the session's own report for that diagnosis; this
script builds nothing and reads no `bin/eqdyna` at all, so it cannot be
fooled by a stale binary the way the by-hand acceptance run was).

Run: python3 testsys/parity/test_standalone_acceptance.py
     (wired into `python3 testsys/run.py accept`, opt-in like parity/perf
     -- needs scripts/create.newcase + case.setup to run cleanly, and the
     committed test.reference.results/ trees, neither guaranteed present
     in every checkout state).
"""
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

TESTSYS = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(TESTSYS))
THRESHOLD = 1e-3  # PROJECT_RULES rule 5's one calibrated tolerance -- never a new number here.

CASES = ('test.tpv8', 'test.tpv104', 'test.tpv1053d')


def sh(cmd, cwd=None, env=None):
    print('+', ' '.join(cmd), (f'(cwd={cwd})' if cwd else ''))
    r = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-4000:]); print(r.stderr[-4000:])
        raise SystemExit(f'FAIL: {" ".join(cmd)} exited {r.returncode}')
    return r


def make_serial_case(case_name, case_dir):
    """create.newcase + force nx=ny=nz=1 + case.setup. No Fortran binary
    involved -- case.setup only writes the bFile/netCDF case-input files
    the standalone solver reads natively."""
    env = dict(os.environ, EQDYNAROOT=REPO_ROOT,
               PYTHONPATH=os.path.join(REPO_ROOT, 'scripts'))
    sh([sys.executable, os.path.join(REPO_ROOT, 'scripts', 'create.newcase'),
        case_dir, case_name], env=env)
    params_path = os.path.join(case_dir, 'user_defined_params.py')
    with open(params_path) as f:
        s = f.read().rstrip('\n')
    with open(params_path, 'w') as f:
        f.write(s + '\n\npar.nx = 1\npar.ny = 1\npar.nz = 1\n')
    sh([sys.executable, 'case.setup'], cwd=case_dir, env=env)


def run_standalone(case_dir):
    """python3 -m eqdyna.standalone <case_dir>, default (jax) backend --
    the exact command a real user would type. Returns frt.txt0's path."""
    env = dict(os.environ, EQDYNAROOT=REPO_ROOT,
               PYTHONPATH=os.path.join(REPO_ROOT, 'python'))
    sh([sys.executable, '-m', 'eqdyna.standalone', case_dir], cwd=REPO_ROOT, env=env)
    frt_path = os.path.join(case_dir, 'frt.txt0')
    if not os.path.isfile(frt_path):
        raise SystemExit(f'FAIL: {frt_path} was not written by eqdyna.standalone')
    return frt_path


def coordinate_aligned_diff(case_name, py_frt_path):
    """dedupe test.reference.results/<case>'s frt.txt0+frt.txt2 by rounded
    coords, lexsort both against the python output by (x,y,z), return
    (max_abs_diff, ok). Raises loudly if the reference tree/files are
    missing -- never a silent skip (PROJECT_RULES rule 2)."""
    ref_dir = os.path.join(REPO_ROOT, 'test.reference.results', case_name)
    ref0_path = os.path.join(ref_dir, 'frt.txt0')
    ref2_path = os.path.join(ref_dir, 'frt.txt2')
    for p in (ref0_path, ref2_path):
        if not os.path.isfile(p):
            raise FileNotFoundError(f'missing reference file {p}')

    ref = np.vstack([np.loadtxt(ref0_path), np.loadtxt(ref2_path)])
    key = np.round(ref[:, :3], 6)
    _, idx = np.unique(key, axis=0, return_index=True)
    ref_deduped = ref[np.sort(idx)]

    py = np.loadtxt(py_frt_path)

    def lexsort_by_xyz(a):
        return a[np.lexsort((a[:, 2], a[:, 1], a[:, 0]))]

    ref_s = lexsort_by_xyz(ref_deduped)
    py_s = lexsort_by_xyz(py)

    if ref_s.shape != py_s.shape:
        raise AssertionError(f'{case_name}: shape mismatch after dedupe/lexsort -- '
                              f'reference {ref_s.shape} vs python {py_s.shape}')
    coord_diff = np.max(np.abs(ref_s[:, :3] - py_s[:, :3]))
    if coord_diff != 0.0:
        raise AssertionError(f'{case_name}: node-coordinate mismatch after lexsort '
                              f'(max abs diff {coord_diff:e}) -- alignment failed, '
                              f'this is not a numeric-tolerance issue')

    max_abs_diff = np.max(np.abs(ref_s - py_s))
    ok = bool(np.allclose(ref_s, py_s, atol=THRESHOLD, rtol=THRESHOLD))
    return max_abs_diff, ok


def run_case(case_name):
    tmp = tempfile.mkdtemp(prefix='eqdyna_accept_' + case_name.replace('.', '_') + '_')
    try:
        case_dir = os.path.join(tmp, case_name)
        make_serial_case(case_name, case_dir)
        frt_path = run_standalone(case_dir)
        max_abs_diff, ok = coordinate_aligned_diff(case_name, frt_path)
        label = 'SUCCESS' if ok else 'FAIL'
        print(f'{label} {case_name} max_abs_diff={max_abs_diff:e} threshold={THRESHOLD:e}')
        return ok
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    print('accept: running %d case(s) serially, one at a time' % len(CASES))
    results = {}
    for case_name in CASES:
        print(f'\n-- {case_name} --')
        try:
            results[case_name] = run_case(case_name)
        except Exception as e:
            print(f'FAIL {case_name} raised: {e!r}')
            results[case_name] = False

    print()
    failures = [c for c, ok in results.items() if not ok]
    if failures:
        print('accept: FAIL -', len(failures), 'case(s) failed:', failures)
        return 1
    print('SUCCESS accept (%d/%d cases)' % (len(results), len(CASES)))
    return 0


if __name__ == '__main__':
    sys.exit(main())

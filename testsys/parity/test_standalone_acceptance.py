#! /usr/bin/env python3
"""
Standalone-Python-EQdyna END-TO-END acceptance tier (`testsys/run.py accept`).

Makes the tpv8/tpv104/tpv1053d acceptance numbers previously recorded ONLY
in a commit message (M7.5/M8 landing commit 057e589; PROJECT_RULES rule 2
"no silent placeholder"/rule 4 "only fresh runs are evidence" -- a number
that only a human can re-derive by reading git log is neither) into a
reproducible, gated, re-runnable check.

For each of test.tpv8 (friclaw=1, planar), test.tpv104 (friclaw=4,
planar), test.tpv1053d (friclaw=5, planar), test.tpv10 (friclaw=1,
insertFaultType=1 dipping planar fault -- Milestone 9), IN SEQUENCE
(never concurrently -- this is a shared machine, same constraint as
testsys/parity/run_parity.py and make_fixtures.py):

  (a) `create.newcase` + `case.setup` a SERIAL (nx=ny=nz=1) case in a
      fresh temp directory (no Fortran binary needed here -- the
      standalone solver reads only the case-input FILES case.setup
      writes, never runs eqdyna).
  (b) `python3 -m eqdyna.standalone <case_dir>` with the DEFAULT (jax)
      backend -- the same invocation a real user would type.
  (c) Coordinate-aligned comparison against the COMMITTED
      test.reference.results/<case>'s frt.txt* files (globbed per case,
      NOT hardcoded to frt.txt0/frt.txt2 -- drv.a6's reference uses
      frt.txt1/frt.txt3; a serial run's frt.txt0 uses a different domain
      decomposition and therefore different node order/count than any
      4-rank reference, so a flat positional diff cannot compare the two
      directly -- see python/eqdyna/standalone/main.py's module
      docstring): dedupe the combined reference by rounded (x,y,z)
      coordinates (partition-boundary nodes are written by more than one
      rank), lexsort both the deduped reference and the python output by
      (x,y,z), then gate on the STRICT ABSOLUTE max diff against a
      per-case bound (PER_CASE_ABS_BOUND below) -- NOT np.allclose (see
      coordinate_aligned_diff's docstring for why: allclose's atol+
      rtol*|b| formula can pass a real scaling bug on large-valued
      columns while the printed max-abs number looks like a failure --
      audit-flagged and fixed here, the printed number and the pass/fail
      decision are now the SAME computed value). THRESHOLD=1e-3 remains
      PROJECT_RULES rule 5's one calibrated OUTER sanity bound; every
      per-case bound is tighter than it, never looser.
  (d) prints `SUCCESS <case> max_abs_diff=<n>` or `FAIL <case> ...` and
      accumulates a non-zero exit if ANY case fails -- no case is allowed
      to fail silently or be skipped without being reported as a failure.

Verified fresh in the session that authored this file: tpv8 8.23e-11,
tpv104 3.39e-7, tpv1053d 6.65e-6 (the last of these AFTER discovering and
fixing a stale `bin/eqdyna` that had pre-dated the fric_tp_h fix -- see
pathway_forward.md/the session's own report for that diagnosis; this
script builds nothing and reads no `bin/eqdyna` at all, so it cannot be
fooled by a stale binary the way the by-hand acceptance run was). All
three numbers are ABSOLUTE max diffs (see PER_CASE_ABS_BOUND).

Run: python3 testsys/parity/test_standalone_acceptance.py
     (wired into `python3 testsys/run.py accept`, opt-in like parity/perf
     -- needs scripts/create.newcase + case.setup to run cleanly, and the
     committed test.reference.results/ trees, neither guaranteed present
     in every checkout state).
"""
import glob
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

TESTSYS = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(TESTSYS))
THRESHOLD = 1e-3  # PROJECT_RULES rule 5's one calibrated tolerance -- never a new number here.

CASES = ('test.tpv8', 'test.tpv104', 'test.tpv1053d', 'test.tpv10')


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


# Per-case ABSOLUTE-diff pass bounds (audit fix -- see module docstring
# note below): each is ~30-100x the observed roundoff for that case (tpv8
# 8.23e-11 -> 1e-8; tpv104 3.39e-7 -> 1e-5; tpv1053d 6.65e-6 -> 1e-4),
# still 5+ orders of magnitude below PROJECT_RULES's outer sanity bound
# (THRESHOLD=1e-3). A case with no recorded observation yet falls back to
# THRESHOLD itself -- never silently tighter or looser than the one
# calibrated number rule 5 allows.
PER_CASE_ABS_BOUND = {
    'test.tpv8': 1e-8,
    'test.tpv104': 1e-5,
    'test.tpv1053d': 1e-4,
    'test.tpv10': 1e-6,  # observed 3.02e-8 (Milestone 9: dipping fault, insertFaultType==1)
}


def coordinate_aligned_diff(case_name, py_frt_path):
    """dedupe test.reference.results/<case>'s reference frt.txt files
    (globbed per-case, NOT hardcoded to frt.txt0/2 -- e.g. drv.a6's
    reference uses frt.txt1/frt.txt3) by rounded coords, lexsort both
    against the python output by (x,y,z). Returns (max_abs_diff, ok, ok
    is a STRICT ABSOLUTE-diff check against PER_CASE_ABS_BOUND -- NOT
    np.allclose (audit fix: allclose's atol+rtol*|b| formula can pass a
    demonstrated scaling bug on large-valued columns, e.g. Tn/Ts/Td ~1e8
    Pa, where rtol*|b| alone is ~1e5 -- decoupling the printed max-abs
    number from the actual pass/fail gate is exactly the silent-fallback
    pattern PROJECT_RULES rule 2 forbids). Raises loudly if the reference
    tree/files are missing -- never a silent skip."""
    ref_dir = os.path.join(REPO_ROOT, 'test.reference.results', case_name)
    ref_files = sorted(glob.glob(os.path.join(ref_dir, 'frt.txt*')))
    if not ref_files:
        raise FileNotFoundError(f'no frt.txt* reference files found under {ref_dir}')

    ref = np.vstack([np.loadtxt(p) for p in ref_files])
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
                              f'reference {ref_s.shape} vs python {py_s.shape} '
                              f'(reference files used: {[os.path.basename(p) for p in ref_files]})')
    # COORD_TOL: for a planar fault (tpv8/104/1053d), node coordinates are
    # exact multiples of the grid spacing on both sides and match to 0.0
    # exactly. For a dipping/rough fault (tpv10/drv.a6), coordinates are the
    # OUTPUT of a real floating-point computation (insertFaultInterface's
    # y-blend) done independently by Fortran and this port -- a few ULPs of
    # difference (~1e-12, the same floor M1's meshCoor check already uses,
    # 1e-9) is expected roundoff, not a misalignment. A genuine node-to-
    # node misalignment would show a diff on the order of the actual grid
    # spacing (meters), many orders of magnitude above this floor -- so a
    # small, fixed tolerance here cannot mask a real alignment bug.
    coord_diff = np.max(np.abs(ref_s[:, :3] - py_s[:, :3]))
    if coord_diff > 1e-9:
        raise AssertionError(f'{case_name}: node-coordinate mismatch after lexsort '
                              f'(max abs diff {coord_diff:e}) -- alignment failed, '
                              f'this is not a numeric-tolerance issue')

    abs_diff = np.abs(ref_s - py_s)
    max_abs_diff = np.max(abs_diff)
    zero_cols = int(np.sum(np.all(np.abs(ref_s) == 0.0, axis=0)))
    print(f'{case_name}: {zero_cols} of {ref_s.shape[1]} reference columns are identically zero '
          f'(unused-for-this-friclaw columns, e.g. fric(20)/fric(23) for friclaw==1)')

    bound = PER_CASE_ABS_BOUND.get(case_name, THRESHOLD)
    ok = bool(max_abs_diff <= bound)
    return max_abs_diff, ok


def run_case(case_name):
    tmp = tempfile.mkdtemp(prefix='eqdyna_accept_' + case_name.replace('.', '_') + '_')
    try:
        case_dir = os.path.join(tmp, case_name)
        make_serial_case(case_name, case_dir)
        frt_path = run_standalone(case_dir)
        max_abs_diff, ok = coordinate_aligned_diff(case_name, frt_path)
        bound = PER_CASE_ABS_BOUND.get(case_name, THRESHOLD)
        label = 'SUCCESS' if ok else 'FAIL'
        print(f'{label} {case_name} max_abs_diff={max_abs_diff:e} bound={bound:e} '
              f'(outer sanity threshold={THRESHOLD:e})')
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

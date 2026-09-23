#! /usr/bin/env python3
"""
The ONE comparison in this repo. Every backend's output is compared here, the
same way, against the same reference, at the same per-case bound.

There used to be two: check.test.py compared a Fortran run's per-rank
frt.txt0..3 positionally against per-rank reference files, and
test_standalone_acceptance.py compared the Python standalone's single
whole-domain frt.txt0 against a deduped, lexsorted concatenation of those same
reference files. Two implementations of "does this match?" is two things to
keep honest, and they did not agree on what a result even was: one of them
carried npx*npy*npz in its file list.

Both are this module now. The comparison is:

    run output  --canonicalise-->  (one row per fault node, ordered by x,y,z)
    reference   (already canonical, committed as frt.canonical.txt)
    ----------------------------------------------------------------
    aligned elementwise -> the CASE's gate (testsys/matrix.py)

Canonicalisation (testsys/frt_canonical.py) is what makes a result a statement
about the physics instead of about the decomposition, so a 4-rank Fortran run,
a serial Fortran run and a serial Python run are all comparable to each other
and to one committed artifact.

TWO GATES, because there are two kinds of case, not because there are two
kinds of backend:
  abs-max      -- strict absolute max |diff|, NOT np.allclose. allclose's
                  atol + rtol*|b| can pass a real scaling bug on the ~1e8 Pa
                  stress columns while the printed max-abs number looks like a
                  failure; decoupling the printed number from the pass/fail
                  decision is the silent-fallback pattern rule 2 forbids. Here
                  the printed number IS the gated number.
  flip-budget  -- test.drv.a6 only: bulk agreement plus an explicit rupture-
                  flip budget, because that case's rupture arrivals are
                  genuinely bistable (see matrix.DRV_A6).

Nothing in here skips. A missing reference, a missing run artifact, a node-set
mismatch and an out-of-bound diff are all failures, and each says which it was.
"""
import os

import numpy as np

from testsys import frt_canonical, matrix

REPO_ROOT = matrix.REPO_ROOT
REFERENCE_ROOT = os.path.join(REPO_ROOT, 'test.reference.results')
CANONICAL_NAME = 'frt.canonical.txt'
NC_NAME = 'fault.dyna.r.nc'

# The gate-term reference, for a case whose own full par.term is not already
# matrix.GATE_TERM_S (see canonical_reference_name below). Chosen name, not
# generated here (rule 7: each reference is its own reviewed commit) -- the
# owner generates it; until then, the cases that need it fail closed.
GATE_TERM_NAME = 'frt.canonical.term5.txt'


def canonical_reference_name(case, term='full'):
    """Which frt reference file name a cell's TERM axis selects.

    term='full' (the default -- unchanged behaviour from before this axis
    existed) always names the one full-length reference, CANONICAL_NAME.

    term='gate' names CANONICAL_NAME too, IF this case's own committed
    par.term (matrix.CASE_FULL_TERM_S[case]) already equals matrix.GATE_TERM_S
    -- nothing is lost by reusing the same file. Otherwise it names
    GATE_TERM_NAME, a second reference this function does NOT create.

    Fail-closed lives downstream, in reference_path(): a case that needs
    GATE_TERM_NAME and does not have it yet raises FileNotFoundError naming
    the missing file. This function never falls back to CANONICAL_NAME for a
    case whose gate term differs from its full term -- that silent
    mismatch (comparing a 5 s run against a 20 s reference) is exactly what
    the owner-approved design forbids.
    """
    if term not in ('gate', 'full'):
        raise ValueError('unknown term %r (expected "gate" or "full")' % term)
    if term == 'full':
        return CANONICAL_NAME
    if matrix.CASE_FULL_TERM_S[case] == matrix.GATE_TERM_S:
        return CANONICAL_NAME
    return GATE_TERM_NAME


# --------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------
def reference_path(case, name=CANONICAL_NAME):
    """The committed reference artifact for a case. Raises if absent: a
    missing reference is a broken checkout, never an excuse to pass."""
    p = os.path.join(REFERENCE_ROOT, case, name)
    if not os.path.isfile(p):
        raise FileNotFoundError(
            'no reference %s for %s at %s -- the reference tree is incomplete; '
            'this is a failure, not a case to skip' % (name, case, p))
    return p


def load_reference(case, term='full'):
    """The case's canonical reference array, (nftnd, 22), at the requested
    TERM axis (default 'full' -- unchanged from before the axis existed).

    ONE file per (case, term), with no rank count in its name. That is the
    point: the gate never globs frt.txt*, so it never encodes how many ranks
    produced the reference. Regenerate with
    `python3 -m testsys.frt_canonical <case_dir>` and commit it deliberately
    (rule 7)."""
    arr = np.loadtxt(reference_path(case, canonical_reference_name(case, term)))
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    if arr.shape[1] != frt_canonical.FRT_COLUMNS:
        raise ValueError('reference %s has %d columns, expected %d'
                         % (reference_path(case, canonical_reference_name(case, term)),
                            arr.shape[1], frt_canonical.FRT_COLUMNS))
    return arr


def load_run(run_dir):
    """The canonical array for a completed run directory, whatever its rank
    count -- 1 file (serial Python), 2 or 4 (MPI Fortran), all the same."""
    return frt_canonical.canonical_from_case(run_dir)


def load_coordinate_aligned(case, run, term='full'):
    """(reference, run) as row-aligned canonical arrays in the same node order,
    at the requested TERM axis (default 'full').

    `run` is a run directory or a single frt-format file. Raises if the two do
    not describe the same fault node set -- that is a different discretisation,
    not a tolerance question."""
    run_arr = (np.loadtxt(run) if os.path.isfile(run) else load_run(run))
    if run_arr.ndim == 1:
        run_arr = run_arr.reshape(1, -1)
    return frt_canonical.align(load_reference(case, term), run_arr)


def align_two_frt_files(path_a, path_b):
    """Align any two frt-format files, NEITHER of which need be the committed
    reference -- e.g. (fresh serial Fortran, standalone Python). Used by
    testsys/parity/evidence_drv_a6_chaos.py so it reuses this alignment rather
    than duplicating it (rule 1)."""
    return frt_canonical.align(np.loadtxt(path_a), np.loadtxt(path_b))


# --------------------------------------------------------------------------
# gate 1: strict absolute max diff
# --------------------------------------------------------------------------
def abs_max_diff(case, ref_a, run_a):
    """(max_abs_diff, ok, index_of_max) for two aligned canonical arrays.

    The printed number and the pass/fail decision are the SAME computed value,
    by construction -- that is the whole reason this is not np.allclose."""
    diff = np.abs(ref_a - run_a)
    idx = np.unravel_index(int(np.argmax(diff)), diff.shape)
    worst = float(diff[idx])
    return worst, bool(worst <= matrix.CASE_BOUND[case]), idx


def abs_max_gate(case, ref_a, run_a):
    """(ok, lines). Strict absolute max |diff| against the case's one bound.

    Reports WHERE the max is (row, column, both values), so a failure is
    diagnosable from the log without re-running anything."""
    bound = matrix.CASE_BOUND[case]
    worst, ok, idx = abs_max_diff(case, ref_a, run_a)
    lines = ['max|diff|=%.6e bound=%.1e at row %d col %d '
             '(ref %.9e vs run %.9e), %d fault nodes compared'
             % (worst, bound, idx[0], idx[1], ref_a[idx], run_a[idx],
                ref_a.shape[0])]
    if not ok:
        lines.append('FAILED the bound by %.1fx -- a bound is never moved to '
                     'make this pass (rule 5)' % (worst / bound))
    return ok, lines


# --------------------------------------------------------------------------
# gate 2: test.drv.a6's flip budget
# --------------------------------------------------------------------------
def flip_decomposition(ref_a, run_a):
    """Existence-flip + timing-shift + bulk-agreement decomposition of two
    already-aligned canonical arrays.

    Returns a dict with every field needed to print a full report even on
    failure -- never a bare bool. Symmetric in its two arguments except for
    which side is named in n_only_ref / n_only_run."""
    col = matrix.DRV_A6['fnft_col']
    sentinel = matrix.DRV_A6['rupture_sentinel']
    ruptured_ref = ref_a[:, col] < sentinel
    ruptured_run = run_a[:, col] < sentinel
    both = ruptured_ref & ruptured_run
    only_ref = ruptured_ref & ~ruptured_run
    only_run = ruptured_run & ~ruptured_ref
    n_both = int(both.sum())
    if n_both == 0:
        raise AssertionError(
            'flip_decomposition: zero fault nodes ruptured in BOTH runs -- one '
            'of them failed to nucleate at all. That is not a flip-count '
            'question.')

    dfnft = np.abs(ref_a[both, col] - run_a[both, col])
    timing_shift = dfnft > matrix.DRV_A6['timing_shift_s']
    n_existence_flips = int(only_ref.sum()) + int(only_run.sum())
    n_timing_shifts = int(timing_shift.sum())
    total_flips = n_existence_flips + n_timing_shifts
    median_fnft_diff = float(np.median(dfnft))

    # Bulk agreement is measured only on matched-arrival nodes (ruptured in
    # both AND not already counted as a timing-shift flip), over the 18
    # non-coordinate, non-fnft physics columns.
    matched = both.copy()
    matched[both] &= ~timing_shift
    phys_cols = [c for c in range(ref_a.shape[1]) if c not in (0, 1, 2, col)]
    phys = np.abs(ref_a[matched][:, phys_cols] - run_a[matched][:, phys_cols])
    phys_max_diff = float(np.max(phys)) if phys.size else 0.0

    flip_mask = only_ref | only_run
    flip_mask[both] |= timing_shift

    return dict(
        n_nodes=int(ref_a.shape[0]), n_ruptured_both=n_both,
        n_matched_arrival=int(matched.sum()), n_only_ref=int(only_ref.sum()),
        n_only_run=int(only_run.sum()), n_existence_flips=n_existence_flips,
        n_timing_shifts=n_timing_shifts, total_flips=total_flips,
        median_fnft_diff=median_fnft_diff, phys_max_diff=phys_max_diff,
        ok_flips=total_flips <= matrix.DRV_A6['total_flip_bound'],
        ok_median=median_fnft_diff <= matrix.DRV_A6['median_fnft_bound'],
        ok_phys=phys_max_diff <= matrix.DRV_A6['phys_max_bound'],
        flip_mask=flip_mask)


def flip_budget_gate(case, ref_a, run_a):
    """(ok, lines) for a flip-budget case. A single scalar cannot distinguish
    "a few hundred marginal nodes flipped rupture arrival, everything else
    matched" from "everything drifted", so all three numbers are printed and
    all three are gated."""
    d = flip_decomposition(ref_a, run_a)
    ok = bool(d['ok_flips'] and d['ok_median'] and d['ok_phys'])
    lines = [
        'nodes=%d ruptured-in-both=%d (matched-arrival=%d)'
        % (d['n_nodes'], d['n_ruptured_both'], d['n_matched_arrival']),
        'flips=%d/%d [%s] (existence=%d: ref-only %d, run-only %d; '
        'timing-shift>%gs=%d)'
        % (d['total_flips'], matrix.DRV_A6['total_flip_bound'],
           'ok' if d['ok_flips'] else 'FAIL', d['n_existence_flips'],
           d['n_only_ref'], d['n_only_run'], matrix.DRV_A6['timing_shift_s'],
           d['n_timing_shifts']),
        'median|dfnft|=%.4fs/%gs [%s]  phys_max=%.6e/%.1e [%s]'
        % (d['median_fnft_diff'], matrix.DRV_A6['median_fnft_bound'],
           'ok' if d['ok_median'] else 'FAIL', d['phys_max_diff'],
           matrix.DRV_A6['phys_max_bound'], 'ok' if d['ok_phys'] else 'FAIL'),
    ]
    return ok, lines


def drv_a6_gate(run, term='full'):
    """test.drv.a6's gate against the committed reference, for callers that
    have a run path rather than aligned arrays (evidence_drv_a6_chaos.py).
    Returns (ok, diagnostics). test.drv.a6's own full term already equals
    matrix.GATE_TERM_S, so term='gate' and term='full' select the same file
    here -- the parameter exists for callers that pass it through generically."""
    ref_a, run_a = load_coordinate_aligned('test.drv.a6', run, term)
    d = flip_decomposition(ref_a, run_a)
    return bool(d['ok_flips'] and d['ok_median'] and d['ok_phys']), d


# --------------------------------------------------------------------------
# the per-cell entry point
# --------------------------------------------------------------------------
def compare_frt(case, run, term='full'):
    """(ok, lines) for one cell's canonical frt output, at the CASE's gate --
    the same gate for every backend, which is what makes the sweep's columns
    comparable to each other. `term` selects which committed reference this
    compares against (compare.canonical_reference_name); it never changes the
    bound (matrix.CASE_BOUND is one number per case, not per term)."""
    ref_a, run_a = load_coordinate_aligned(case, run, term)
    if matrix.GATE[case] == 'abs-max':
        return abs_max_gate(case, ref_a, run_a)
    return flip_budget_gate(case, ref_a, run_a)


def compare_nc_files(fn1, fn2, threshold=matrix.THRESHOLD):
    """fault.dyna.r.nc comparison: variable set + attributes must match, and
    every variable must agree within the one calibrated threshold (rule 5).

    Bit-exact data equality is deliberately NOT required -- a parallel MPI
    rerun cannot promise it (reduction order varies) -- but a changed variable
    set or changed attributes is a hard failure, not a tolerance question.
    Returns the printed SUCCESS/FAIL string."""
    from netCDF4 import Dataset

    verdict = 'SUCCESS ' + fn1 + ' ' + fn2

    def attrs(obj):
        return {k: obj.getncattr(k) for k in obj.ncattrs()}

    def attrs_equal(a, b):
        def val_eq(x, y):
            try:
                return np.array_equal(x, y, equal_nan=True)
            except TypeError:      # non-numeric attrs (strings, mixed)
                return np.array_equal(x, y)
        return set(a) == set(b) and all(val_eq(a[k], b[k]) for k in a)

    f1 = Dataset(fn1, 'r')
    f2 = Dataset(fn2, 'r')
    try:
        metadata_equal = (
            set(f1.variables) == set(f2.variables)
            and attrs_equal(attrs(f1), attrs(f2))
            and all(attrs_equal(attrs(f1.variables[v]), attrs(f2.variables[v]))
                    for v in f1.variables)
        )
        for var in f1.variables:
            var1 = f1.variables[var]
            var2 = f2.variables[var]
            if var1.dimensions != var2.dimensions:
                verdict = 'FAIL var dim ' + fn1 + ' ' + fn2
            elif not np.allclose(np.asarray(var1[:]), np.asarray(var2[:]),
                                 rtol=threshold, atol=threshold):
                verdict = 'FAIL var numbers ' + fn1 + ' ' + fn2
        if not metadata_equal and verdict.startswith('SUCCESS'):
            verdict = 'FAIL metadata ' + fn1 + ' ' + fn2
    finally:
        f1.close()
        f2.close()
    return verdict


def compare_nc(case, run_dir):
    """(ok, lines) for the fault.dyna.r.nc artifact of one cell."""
    ref = reference_path(case, NC_NAME)
    run = os.path.join(run_dir, NC_NAME)
    if not os.path.isfile(run):
        return False, ['nc: FAIL missing %s -- plotRuptureDynamics did not run '
                       'or did not write it' % run]
    verdict = compare_nc_files(ref, run)
    return verdict.startswith('SUCCESS'), ['nc: %s (threshold=%.0e)'
                                           % (verdict, matrix.THRESHOLD)]


def compare_cell(case, backend, run_dir, term='full'):
    """(ok, lines) for one (case, backend) cell: every artifact that backend
    produces, compared against the committed reference, at the requested
    TERM axis (default 'full', unchanged from before the axis existed).

    Which artifacts is DATA (matrix.ARTIFACTS), so the cell's report states
    what it covered -- 'frt' alone for the python backends, 'frt+nc' for
    fortran -- with ONE term-specific exception, stated rather than silent:
    the 'nc' artifact's own committed reference (fault.dyna.r.nc) is a
    FULL-TERM-ONLY artifact -- there is no gate-term counterpart planned (the
    owner-approved design adds a second reference only for 'frt'). For a case
    whose gate term differs from its full term, comparing nc at term='gate'
    would be a guaranteed, physics-free mismatch (a truncated run against a
    full-length reference) -- so for THOSE cases only, nc is not compared at
    term='gate', and the omission is printed, never swallowed. A case whose
    full term ALREADY equals matrix.GATE_TERM_S (e.g. test.tpv8, the one case
    CI's smoke job runs) is unaffected: its one committed nc reference IS a
    gate-term reference, so nc stays compared there. This is a term-scoping
    DECISION, not a fallback -- it never substitutes the full-term reference
    for a missing gate-term one."""
    if term not in ('gate', 'full'):
        raise ValueError('unknown term %r (expected "gate" or "full")' % term)
    artifacts = matrix.ARTIFACTS[backend]
    skipped_nc = (term == 'gate' and 'nc' in artifacts
                 and matrix.CASE_FULL_TERM_S[case] != matrix.GATE_TERM_S)
    if skipped_nc:
        artifacts = tuple(a for a in artifacts if a != 'nc')
    ok, lines = True, []
    for artifact in artifacts:
        if artifact == 'frt':
            a_ok, a_lines = compare_frt(case, run_dir, term)
        elif artifact == 'nc':
            a_ok, a_lines = compare_nc(case, run_dir)
        else:
            raise ValueError('unknown artifact %r for backend %r -- every '
                             'artifact needs a comparison, none is skipped'
                             % (artifact, backend))
        ok = ok and a_ok
        lines.extend(a_lines)
    if skipped_nc:
        lines.append('nc: not compared at term=gate -- fault.dyna.r.nc has '
                     'no gate-term reference by design (see compare_cell '
                     'docstring), not a coverage gap')
    return ok, lines

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
import glob
import os
import re
import threading

import numpy as np

from testsys import frt_canonical, matrix

REPO_ROOT = matrix.REPO_ROOT
REFERENCE_ROOT = os.path.join(REPO_ROOT, 'test.reference.results')
CANONICAL_NAME = 'frt.canonical.txt'
NC_NAME = 'fault.dyna.r.nc'
_NC_LOCK = threading.Lock()


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


def load_reference(case):
    """The case's canonical reference array, (nftnd, 22) -- the ONE
    committed frt.canonical.txt, at matrix.GATE_TERM_S (there is no second
    term and no second reference file for any case, 2026-09-23).

    ONE file per case, with no rank count in its name. That is the
    point: the gate never globs frt.txt*, so it never encodes how many ranks
    produced the reference. Regenerate with
    `python3 -m testsys.frt_canonical <case_dir>` and commit it deliberately
    (rule 7)."""
    arr = np.loadtxt(reference_path(case))
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    if arr.shape[1] != frt_canonical.FRT_COLUMNS:
        raise ValueError('reference %s has %d columns, expected %d'
                         % (reference_path(case), arr.shape[1],
                            frt_canonical.FRT_COLUMNS))
    return arr


def load_run(run_dir):
    """The canonical array for a completed run directory, whatever its rank
    count -- 1 file (serial Python), 2 or 4 (MPI Fortran), all the same."""
    return frt_canonical.canonical_from_case(run_dir)


def load_coordinate_aligned(case, run):
    """(reference, run) as row-aligned canonical arrays in the same node order.

    `run` is a run directory or a single frt-format file. Raises if the two do
    not describe the same fault node set -- that is a different discretisation,
    not a tolerance question."""
    run_arr = (np.loadtxt(run) if os.path.isfile(run) else load_run(run))
    if run_arr.ndim == 1:
        run_arr = run_arr.reshape(1, -1)
    return frt_canonical.align(load_reference(case), run_arr)


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


def drv_a6_gate(run):
    """test.drv.a6's gate against the committed reference, for callers that
    have a run path rather than aligned arrays (evidence_drv_a6_chaos.py).
    Returns (ok, diagnostics)."""
    ref_a, run_a = load_coordinate_aligned('test.drv.a6', run)
    d = flip_decomposition(ref_a, run_a)
    return bool(d['ok_flips'] and d['ok_median'] and d['ok_phys']), d


# --------------------------------------------------------------------------
# the per-cell entry point
# --------------------------------------------------------------------------
def compare_frt(case, run):
    """(ok, lines) for one cell's canonical frt output, at the CASE's gate --
    the same gate for every backend, which is what makes the sweep's columns
    comparable to each other, against the ONE committed reference (there is
    no term axis to select between)."""
    ref_a, run_a = load_coordinate_aligned(case, run)
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
    # The sweep compares cells from concurrent threads, and since python-jax
    # cells carry 'nc' too, two netCDF4/HDF5 opens could overlap: the HDF5
    # library in use is not thread-safe, and the first everyday sweep of this
    # gate died with SIGSEGV (exit -11) in exactly that window. One lock.
    with _NC_LOCK:
        return _compare_nc_files_locked(fn1, fn2, threshold)


def _compare_nc_files_locked(fn1, fn2, threshold):
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
        # Item 106(1): with zero variables the loop above compares nothing
        # and the verdict stayed SUCCESS -- a pass must mean the check ran
        # (rule 2). Two empty files are "nothing to compare", not "equal".
        if not f1.variables and verdict.startswith('SUCCESS'):
            verdict = 'FAIL no variables compared ' + fn1 + ' ' + fn2
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


# --------------------------------------------------------------------------
# station files: the n-stress sign convention (board row 22a)
# --------------------------------------------------------------------------
ONFAULT_STATION_RE = re.compile(r'^faultst(-?\d+)dp(\d+)\.txt$')


_FORTRAN_E_DROPPED = re.compile(r'^([-+]?\d*\.\d*)([-+]\d{3})$')


def _fortran_float(tok):
    """float() of one Fortran-written value. A plain Ew.d edit (the station
    writers' E16.7) DROPS the 'E' when the exponent needs three digits, e.g.
    '0.1234567-100' for 1.234567e-101; float() would raise on it."""
    m = _FORTRAN_E_DROPPED.match(tok)
    return float(m.group(1) + 'E' + m.group(2)) if m else float(tok)


def read_station_file(path):
    """(field_names, data) for one faultst*/body* station file, as written by
    library_output.f90: `#` header lines, one line of field names starting
    with `t`, then one numeric row per time step. Raises on a file with no
    names line or no data -- an unreadable station is a failure, not a
    station to skip."""
    names, rows = None, []
    with open(path) as f:
        for line in f:
            tok = line.split()
            if not tok or tok[0].startswith('#'):
                continue
            if tok[0] == 't':
                names = tok
                continue
            rows.append([_fortran_float(v) for v in tok])
    if names is None or not rows:
        raise ValueError('%s: no field-name line or no data rows' % path)
    data = np.asarray(rows)
    if data.shape[1] != len(names):
        raise ValueError('%s: %d data columns but %d field names'
                         % (path, data.shape[1], len(names)))
    return names, data


def nstress_sign_gate(case, run_dir):
    """(ok, lines): at every BURIED on-fault station (down-dip field of the
    filename > 0; a surface station's normal stress is legitimately zero), the
    n-stress at the first recorded step must carry the sign an initially
    COMPRESSIVE stress has in the case's SCEC convention
    (matrix.NSTRESS_CONVENTION): negative where "positive means extension",
    positive where "positive means compression". Zero is a failure (the sign
    cannot be read), and so is a run with no buried on-fault station at all
    -- a check that examined nothing must not read as a pass (rule 2)."""
    convention, citation = matrix.NSTRESS_CONVENTION[case]
    expected = -1.0 if convention == 'extension' else 1.0
    buried, unparsed = [], []
    for p in sorted(glob.glob(os.path.join(run_dir, 'faultst*.txt'))):
        m = ONFAULT_STATION_RE.match(os.path.basename(p))
        if not m:
            unparsed.append(os.path.basename(p))
        elif int(m.group(2)) > 0:
            buried.append(p)
    if unparsed:
        return False, ['nsign: FAIL on-fault station file name(s) not of the '
                       'form faultst<sss>dp<ddd>.txt, so their depth cannot '
                       'be read: %s' % ', '.join(unparsed)]
    if not buried:
        return False, ['nsign: FAIL no buried on-fault station file '
                       '(faultst*dp<ddd>.txt, ddd > 0) in %s -- nothing to '
                       'check' % run_dir]
    # Every on-fault station file, surface ones included, must carry data. A
    # buried station's initial traction is nonzero; a surface station can start
    # at zero stress (tpv8/10/36/37) and is kept nonzero here by the waves that
    # reach it within the gate term (<= ~19 km from a t=0 nucleation, well
    # inside 5 s). An all-zero file is the phantom
    # station eqdyna3d.f90's allocInitAfterMeshGen used to make a rank with no
    # matched station write (faultst000dp000.txt, requested by nobody).
    empty = []
    for p in sorted(glob.glob(os.path.join(run_dir, 'faultst*.txt'))):
        names, data = read_station_file(p)
        if not np.any(data[:, 1:]):
            empty.append(os.path.basename(p))
    if empty:
        return False, ['nsign: FAIL on-fault station file(s) with no nonzero '
                       'value in any column: %s -- a phantom station, not a '
                       'real one' % ', '.join(empty)]
    values, bad = [], []
    for p in buried:
        names, data = read_station_file(p)
        v = float(data[0, names.index('n-stress')])
        values.append(v)
        if np.sign(v) != expected:
            bad.append('%s n-stress %.6g MPa' % (os.path.basename(p), v))
    lines = ['nsign: positive means %s (%s); %d buried on-fault stations, '
             'first-step n-stress in [%.6g, %.6g] MPa'
             % (convention, citation, len(buried), min(values), max(values))]
    if bad:
        lines.append('nsign: FAIL %d of %d stations carry the wrong sign for '
                     'an initially compressive stress: %s'
                     % (len(bad), len(buried), '; '.join(bad)))
    return not bad, lines


# --------------------------------------------------------------------------
# station gate: on/off-fault time series, normalized (owner design,
# mission "iris/station-gate", 2026-09-24)
# --------------------------------------------------------------------------
STATION_REFERENCE_DIRNAME = 'stations'


def station_reference_path(case, filename):
    """The committed reference station file for one of matrix.GATE_STATIONS'
    selected files. Raises if absent -- same contract as reference_path
    above: a missing reference is a broken checkout, not a skip."""
    p = os.path.join(REFERENCE_ROOT, case, STATION_REFERENCE_DIRNAME, filename)
    if not os.path.isfile(p):
        raise FileNotFoundError(
            'no station reference %s for %s at %s -- the reference tree is '
            'incomplete; this is a failure, not a case to skip' % (filename, case, p))
    return p


def _load_station_pair(case, kind, filename, run_dir):
    """(names, ref_data, run_data) for one selected station file, having
    already checked field names, column count and the time axis match.
    Raises (never returns a partial/interpolated result) on any of: the run
    not writing the file, a field-name/column-count mismatch, or a time-axis
    mismatch (different row count or different t values) -- rule 2, and the
    mission's explicit "never an interpolation"."""
    ref_path = station_reference_path(case, filename)
    run_path = os.path.join(run_dir, filename)
    if not os.path.isfile(run_path):
        raise FileNotFoundError(
            'station gate: %s x %s: %s is a SELECTED station '
            '(matrix.GATE_STATIONS) but the run wrote no such file -- a '
            'selected file missing from a run is a FAIL, not a skip'
            % (case, kind, filename))
    ref_names, ref_data = read_station_file(ref_path)
    run_names, run_data = read_station_file(run_path)
    if ref_names != run_names:
        raise ValueError(
            'station gate: %s %s: field-name line differs -- ref %r vs run '
            '%r (%s)' % (case, filename, ref_names, run_names, ref_path))
    if ref_data.shape[0] != run_data.shape[0]:
        raise ValueError(
            'station gate: %s %s: %d reference rows vs %d run rows -- a row '
            'count mismatch is a FAIL, never an interpolation'
            % (case, filename, ref_data.shape[0], run_data.shape[0]))
    if not np.allclose(ref_data[:, 0], run_data[:, 0], atol=1e-9, rtol=0.0):
        bad = int(np.argmax(np.abs(ref_data[:, 0] - run_data[:, 0])))
        raise ValueError(
            'station gate: %s %s: time column mismatch at row %d (ref t=%.6e '
            's vs run t=%.6e s) -- a FAIL, never an interpolation'
            % (case, filename, bad, ref_data[bad, 0], run_data[bad, 0]))
    return ref_names, ref_data, run_data


def normalize_station_error(ref_datas, run_datas, col):
    """e_q = max_t|run-ref| / max(S_q, matrix.STATION_ZERO_FLOOR), where
    S_q = max over the given stations (one case, one kind: 'on' or 'off') of
    max_t|ref_q|. `ref_datas`/`run_datas` are lists of per-station (nsteps,
    ncols) arrays, already time-aligned by `_load_station_pair`. Returns
    (e, worst_diff, S_q, station_index_of_worst).

    The max(S_q, FLOOR) clamp is the mission's "S_q == 0" rule generalised to
    a MEASURED floor rather than exact zero -- see matrix.STATION_ZERO_FLOOR's
    own docstring for the measurement that motivated it (a genuinely-zero
    physical component never lands at literal 0.0 in a floating-point run;
    it lands at that run's own roundoff noise, and two roundoff numbers
    divided by each other is not a physics statement)."""
    S_q = max(float(np.max(np.abs(d[:, col]))) for d in ref_datas)
    diffs = [float(np.max(np.abs(r[:, col] - u[:, col])))
             for r, u in zip(ref_datas, run_datas)]
    worst_i = int(np.argmax(diffs))
    worst_diff = diffs[worst_i]
    scale = max(S_q, matrix.STATION_ZERO_FLOOR)
    return worst_diff / scale, worst_diff, S_q, worst_i


def station_gate(case, run_dir):
    """(ok, lines) for one cell's station artifact.

    Reads every file matrix.GATE_STATIONS[case] selects (both 'on' and
    'off'), checks field names/column count/time axis (rule 2: never an
    interpolation), computes the case-level scale S_q per (kind, column) and
    the normalized error e_q = diff/max(S_q, FLOOR), and gates every e_q
    against the ONE matrix.STATION_BOUND[case] (rule 5). A cell declared in
    matrix.STATION_UNSUPPORTED is handled by compare_cell, which never calls
    this for it."""
    bound = matrix.STATION_BOUND[case]
    kinds = matrix.GATE_STATIONS[case]
    worst = (0.0, None)
    lines, bad_cols = [], []
    for kind, files in kinds.items():
        names = None
        ref_datas, run_datas = [], []
        for fn in files:
            try:
                n, ref_d, run_d = _load_station_pair(case, kind, fn, run_dir)
            except (OSError, ValueError) as e:
                # A missing selected file, a field-name/column mismatch or a
                # time-axis mismatch is a FAILED cell with its reason, the
                # same way an out-of-bound diff is -- never a crash that
                # reads differently in the sweep summary.
                return False, ['station: FAIL %s' % e]
            names = n
            ref_datas.append(ref_d)
            run_datas.append(run_d)
        ncols = ref_datas[0].shape[1]
        for col in range(1, ncols):
            e, worst_diff, S_q, worst_i = normalize_station_error(
                ref_datas, run_datas, col)
            floor_used = S_q < matrix.STATION_ZERO_FLOOR
            lines.append(
                '  %-3s %-16s S_q=%.4e%s diff=%.4e e=%.4e/%.1e %s'
                % (kind, names[col], S_q,
                   ' (FLOOR-clamped)' if floor_used else '', worst_diff, e,
                   bound, 'ok' if e <= bound else 'FAIL'))
            if not e <= bound:
                # `not e <= bound`, never `e > bound`: a NaN anywhere in a
                # run's column makes e NaN, and NaN > bound is False -- the
                # form that read a NaN station as green.
                bad_cols.append('%s %s at %s' % (kind, names[col], files[worst_i]))
            if e > worst[0] or np.isnan(e):
                worst = (e, '%s %s at %s' % (kind, names[col], files[worst_i]))
    ok = not bad_cols
    header = ('station: %d on-fault + %d off-fault file(s), worst e=%.4e '
              'bound=%.1e (%s)%s'
             % (len(kinds['on']), len(kinds['off']), worst[0], bound,
                worst[1], '' if ok else ' -- FAILED'))
    return ok, [header] + lines


def compare_cell(case, backend, run_dir):
    """(ok, lines) for one (case, backend) cell: every artifact that backend
    produces, compared against the ONE committed reference (there is no term
    axis and no second reference file, 2026-09-23).

    Which artifacts is DATA (matrix.ARTIFACTS), so the cell's report states
    what it covered -- 'frt' alone for the python backends, 'frt+nc+nsign'
    for fortran. Every declared artifact is compared, unconditionally; nothing
    here is ever skipped."""
    artifacts = matrix.ARTIFACTS[backend]
    ok, lines = True, []
    for artifact in artifacts:
        if artifact == 'frt':
            a_ok, a_lines = compare_frt(case, run_dir)
        elif artifact == 'nc':
            if (case, backend) in matrix.NC_UNSUPPORTED:
                a_ok, a_lines = True, ['nc: DECLARED UNSUPPORTED for %s x %s -- %s'
                                       % (case, backend,
                                          matrix.NC_UNSUPPORTED[(case, backend)])]
            else:
                a_ok, a_lines = compare_nc(case, run_dir)
        elif artifact == 'nsign':
            a_ok, a_lines = nstress_sign_gate(case, run_dir)
        elif artifact == 'station':
            # A sub-cell declaration (matrix.STATION_UNSUPPORTED, keyed per
            # (case, backend)): reported, not silently absent, and it does NOT
            # fail the cell -- test.drv.a6 x python-jax still passes or fails
            # on its frt/nsign gates. Not the same state machine as
            # matrix.UNSUPPORTED (whole-CELL unsupported): the cell runs and
            # is gated on everything else it carries; only this artifact, for
            # this one cell, is declared, with the measured reason printed.
            if (case, backend) in matrix.STATION_UNSUPPORTED:
                a_ok, a_lines = True, [
                    'station: DECLARED UNSUPPORTED for %s x %s -- %s'
                    % (case, backend, matrix.STATION_UNSUPPORTED[(case, backend)])]
            else:
                a_ok, a_lines = station_gate(case, run_dir)
        else:
            raise ValueError('unknown artifact %r for backend %r -- every '
                             'artifact needs a comparison, none is skipped'
                             % (artifact, backend))
        ok = ok and a_ok
        lines.extend(a_lines)
    return ok, lines

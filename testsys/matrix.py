#! /usr/bin/env python3
"""
The e2e sweep's TABLE: which (case, backend) cells exist, which are supported,
and the ONE bound each case is gated at.

There is one test in this repo -- the e2e sweep -- and it is two nested loops:

    for case in CASES:
        for backend in BACKENDS:
            run(case, backend) -> canonical frt -> compare against the ONE
                                  committed reference for that case, at THAT
                                  CASE's bound

`backend` is an axis of the sweep, not a tier. What used to be the separate
`accept` tier (the same standalone solver, the same reference tree, a second
alignment+comparison implementation, and a case list of 5 while e2e gated 8 --
so "SUCCESS accept (5/5 cases)" meant five of five LISTED, jax only, i.e. 5 of
16 case x backend combinations) is now two columns of this table. Adding a GPU
backend is appending to BACKENDS, not writing a new runner.

WHY THE TABLE IS DATA, NOT CONTROL FLOW
A cell is either SUPPORTED -- it runs, and it is gated -- or it is DECLARED
UNSUPPORTED here with a recorded reason. There is no third state, no skip, and
no silent absence. PROJECT_RULES rule 2 requires "I could not check this" and
"this is fine" to produce different exit codes, so asking the sweep for a cell
this table declares unsupported is a hard FAILURE (see cells()).

WHY THE BOUND BELONGS TO THE CASE, NOT THE BACKEND
A per-backend bound makes the backends incomparable: "jax passed and numpy
failed" would not be a statement about the solver, it would be a statement
about which number each column was allowed. One bound per case, applied
identically to every backend, is what makes a row of this sweep quantitative.

The bounds below are UNCHANGED from the accept tier they came from (rule 5:
one calibrated definition of pass; and no bound is ever moved to make a cell
go green). Each is 15-120x the roundoff observed for that case on the jax
backend, and every one is at or below THRESHOLD, the repo's single outer
sanity bound.
"""
import os
import sys

TESTSYS = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(TESTSYS)
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from testNameList import nameList as _NAME_LIST, coreNumList as _CORE_NUM_LIST

# CASES WITHDRAWN FROM THE SWEEP, kept as reference.
#
# A case listed here is NOT gated: the sweep does not run it, and a green sweep
# does not speak for it. It is not deleted either -- its committed reference
# stays, the unit tier keeps exercising it, and the open question is written
# down where work gets picked up (pathway_forward.md). "Withdrawn with a
# recorded reason and a live reference" is a third declared state, not a skip:
# the coverage report prints it on every run, and naming it explicitly is a
# hard failure (see cells()), so it cannot quietly become forgotten coverage.
REFERENCE_ONLY = {
    'test.drv.a6':
        'Withdrawn 2026-09-15 pending a numerics answer, not a gate answer. '
        'Fortran compared against its OWN 4-rank reference, changing NOTHING '
        'but the decomposition, gives 329 arrival flips of 5151 and a 2.12e7 '
        'max difference. A case whose Fortran-vs-Fortran floor is that high '
        'cannot say anything about a backend until we know whether that is '
        'genuine bistability or dx=500 under-resolution with friclaw=4. Its '
        'reference, its flip-budget gate (DRV_A6 below) and the evidence '
        'script all stay; see pathway_forward.md "drv.a6 decomposition '
        'sensitivity" for the two experiments that decide it.',
}

# The case list is testNameList.py's, not a second copy of it (rule 1). A case
# added there with no entry below is a hard import-time failure, not a silently
# ungated case -- see the consistency block at the bottom of this module.
ALL_CASES = tuple(_NAME_LIST)               # every case the repo has a case for
CASES = tuple(c for c in _NAME_LIST if c not in REFERENCE_ONLY)   # gated ones
FORTRAN_RANKS = dict(zip(_NAME_LIST, _CORE_NUM_LIST))

BACKENDS = ('fortran', 'python-numpy', 'python-jax')

# Which artifacts each backend produces, and therefore what gets compared.
#   frt -- the canonical fault-node table (testsys/frt_canonical.py). Every
#          backend produces it; it is what makes the sweep ONE comparison.
#   nc  -- fault.dyna.r.nc, written by scripts/plotRuptureDynamics from the frt
#          output. DELIBERATELY still compared, for the fortran column only:
#          it is strictly lossy relative to frt (12 resampled dip x strike
#          variables vs 18 physics columns at native node resolution), so it
#          adds no physics coverage -- but it is the ONLY thing that exercises
#          plotRuptureDynamics, post-processing that would otherwise be gated
#          by nothing. It is backend DATA here, so a green run states that the
#          python columns did not cover it.
ARTIFACTS = {
    'fortran': ('frt', 'nc'),
    'python-numpy': ('frt',),
    'python-jax': ('frt',),
}

THRESHOLD = 1e-3  # PROJECT_RULES rule 5's one outer sanity bound.

# One bound per case. None means "this case is not gated on a scalar" -- see
# GATE/DRV_A6 below. Every case in CASES must appear here.
CASE_BOUND = {
    'test.tpv8': 1e-8,        # observed 8.23e-11 (python-jax vs reference)
    'test.tpv104': 1e-5,      # observed 3.39e-07
    'test.tpv1053d': 1e-4,    # observed 6.65e-06
    'test.tpv10': 1e-6,       # observed 3.02e-08 (dipping fault, insertFaultType==1)
    'test.meng2023a': THRESHOLD,   # fortran-only column so far; no python run has
    'test.meng2023cb': THRESHOLD,  # ever been measured for these, so they carry
    'test.tpv29': THRESHOLD,       # the outer bound EXPLICITLY rather than a
                                   # tighter number nobody has observed.
}

GATE = {
    'test.tpv8': 'abs-max',
    'test.tpv104': 'abs-max',
    'test.tpv1053d': 'abs-max',
    'test.tpv10': 'abs-max',
    'test.meng2023a': 'abs-max',
    'test.meng2023cb': 'abs-max',
    'test.tpv29': 'abs-max',
    # test.drv.a6 was gated 'flip-budget' here until 2026-09-15; it is now in
    # REFERENCE_ONLY and therefore has no entry. The gate itself is kept below
    # and still implemented in compare.py -- withdrawing the case is not a
    # verdict that the gate was wrong, and whichever way the numerics question
    # resolves, this is what it comes back as.
}

# test.drv.a6's flip-budget gate. NOT WIRED INTO THE SWEEP -- see
# REFERENCE_ONLY. Still live for testsys/parity/evidence_drv_a6_chaos.py, which
# is how these numbers are regenerated; they are measured, not chosen, and must
# not be changed without re-running it fresh (rule 4).
#
# A scalar max-abs bound cannot distinguish "a few hundred marginal nodes
# flipped" from "everything drifted a little", which is why this case was given
# bulk agreement PLUS an explicit flip budget rather than a looser threshold.
DRV_A6 = {
    'fnft_col': 3,               # canonical frt column layout: x, y, z, fnft, ...
    'rupture_sentinel': 1.0e4,   # fnft sentinel is 99999.0 in Fortran
                                 # (eqdyna3d.f90:139) and in the port; below
                                 # 1e4 means "this node ruptured".
    'timing_shift_s': 1.0,       # a ruptured-in-both node counts as a flip only
                                 # past this |delta fnft|.
    'median_fnft_bound': 0.1,    # seconds; measured 0.0417 s in all three of
                                 # A/B/C -- 10x margin.
    'phys_max_bound': 1.1e9,     # 30x the measured 3.996e7 ceiling over the 18
                                 # non-coordinate/non-fnft columns, restricted
                                 # to matched-arrival nodes (those columns carry
                                 # ~1e8 Pa stress fields).
    'total_flip_bound': 450,     # A=372 (Fortran-only decomposition floor),
                                 # B=439 (python vs serial fortran). 450 is
                                 # minimal headroom above the larger measured
                                 # number, fixed BEFORE C (391) was measured.
                                 # A LATER, independent Fortran-vs-Fortran run
                                 # (serial vs the committed 4-rank reference,
                                 # 2026-09-15) measured 329 flips, max diff
                                 # 2.12e7, p50 4.03e-03, p99 5.04e+05, and a
                                 # median arrival difference of 0.0417 s =
                                 # exactly one time step. Recorded, not
                                 # substituted: A and this run are different
                                 # experiments and 450 stays where it was set.
                                 # That 329 is the floor imposed by
                                 # decomposition ALONE is what withdrew the
                                 # case -- python-jax 397 and python-numpy 461
                                 # straddle a budget whose own floor is 329.
}

# DECLARED UNSUPPORTED CELLS. A reason, verified in the source, not a guess.
# Absence from this table means "supported"; presence means the sweep refuses
# to pretend it covered the cell.
UNSUPPORTED = {
    # EMPTY, and that is a measured statement, not an oversight.
    #
    # It held four cells -- test.meng2023a and test.meng2023cb on both python
    # backends -- with the reason "the Python port has no solveSWTW path".
    # That reason is now FALSE: eqdyna/faulting.py has solveSWTW, dispatched at
    # `if friclaw <= 2` exactly where faulting.f90:17 dispatches it, and
    # eqdyna/fric.py has time_weak (slip_weak with slip/SW_D0 replaced by
    # trupt/TW_T0; trupt = timeElapsed - fnft, and the 99999.0 sentinel makes
    # trupt negative for an unruptured node so the fs arm fires with no special
    # case). Leaving a stale entry here would be a lie about coverage.
    #
    # First python runs these two cases have ever had, against their 1e-3
    # bound -- NOT tightened here; a bound is calibrated in its own commit,
    # never in the change that first makes the cell green:
    #     meng2023a  x python-numpy  4.628301e-09
    #     meng2023a  x python-jax    3.150106e-09
    #     meng2023cb x python-numpy  5.054473e-09
    #     meng2023cb x python-jax    5.042553e-09
    # Both cases also set par.tpv = 201, so they exercise swtwNucleation's
    # smoothed forced-rupture branch at the same time.
}

# Cells CI runs, and the measured reason the rest are left out. A GitHub
# ubuntu-22.04 runner has 7 GB and is already holding the Fortran build;
# exceeding it produced a SIGTERM with no output (exit 143), which is a
# resource kill and not a test result. This is a declared, measured decision
# with the numbers next to it -- not an env var whose value nobody can audit.
MEASURED_PEAK_RSS_GB = {
    # Full-length runs, one case at a time on the 64-core development box.
    # The two CI cells were measured with `/usr/bin/time -v` (exact peak):
    ('test.tpv8', 'python-numpy'): 1.45,   # 62.5 s wall
    ('test.tpv8', 'python-jax'): 1.91,     # 30.7 s wall
    # The rest are from earlier per-case measurements on the same box, kept
    # because they are what the CI exclusion rests on:
    ('test.tpv104', 'python-jax'): 3.42,
    ('test.tpv10', 'python-jax'): 3.23,
    ('test.tpv1053d', 'python-jax'): 3.95,
    # drv.a6 x python-jax was re-measured during the full sweep by sampling
    # the process RSS every 20 s: 9.57 GB, i.e. a LOWER BOUND on its true
    # peak and well above the earlier 7.40 GB figure. Either way this cell
    # cannot run on a 7 GB runner; the newer, larger number is recorded
    # because that is the one a widening decision must be made against.
    ('test.drv.a6', 'python-jax'): 9.57,
}
CI_RUNNER_RAM_GB = 7.0
CI_CELLS = (
    tuple((c, 'fortran') for c in CASES)
    + (('test.tpv8', 'python-numpy'), ('test.tpv8', 'python-jax'))
)
# WHAT THIS LIST LEAVES OUT, AND WHY -- said here rather than implied:
#   * python cells for tpv10/tpv104/tpv1053d: measured 3.2-4.0 GB above,
#     against a 7 GB runner already holding the Fortran build. A resource kill.
#   * python cells for tpv29: these FAIL today (both backends, identical
#     max|diff|=7.10e7 -- see the measured sweep). They are not omitted for
#     memory, they are omitted because CI would be permanently red on a known
#     solver gap, and a permanently red gate stops being read. The full sweep
#     (`testsys/run.py e2e`) runs them and fails, which is where that failure
#     is supposed to be visible. Moving a bound to make them pass is not an
#     option (rule 5).
#   * every drv.a6 cell, including the fortran one: the case is in
#     REFERENCE_ONLY, so it is out of CASES and out of this list with it.
# A green CI run therefore means 9 of 21 cells, and says so.


def is_supported(case, backend):
    return (case, backend) not in UNSUPPORTED


def is_reference_only(case):
    """True for a case kept as reference but withdrawn from the sweep."""
    return case in REFERENCE_ONLY


def reference_only_reason(case):
    """Why a case was withdrawn. Raises for a gated case: asking why a gated
    case was skipped is itself a bug (rule 2)."""
    if case not in REFERENCE_ONLY:
        raise KeyError('%s is gated -- it has no reference-only reason' % case)
    return REFERENCE_ONLY[case]


def unsupported_reason(case, backend):
    """The recorded reason a cell is not gated. Raises for a supported cell:
    asking why a supported cell was skipped is itself a bug (rule 2)."""
    if is_supported(case, backend):
        raise KeyError('%s x %s is supported -- it has no unsupported reason'
                       % (case, backend))
    return UNSUPPORTED[(case, backend)]


def gate_description(case):
    if GATE[case] == 'abs-max':
        return 'max|diff| <= %.1e' % CASE_BOUND[case]
    return ('flips <= %d, median|dfnft| <= %.2gs, phys <= %.1e'
            % (DRV_A6['total_flip_bound'], DRV_A6['median_fnft_bound'],
               DRV_A6['phys_max_bound']))


def cells(cases=None, backends=None):
    """Return (runnable, declared_unsupported) for the requested selection.

    `runnable` is the list of (case, backend) the sweep will execute;
    `declared_unsupported` is the list of (case, backend, reason) INSIDE the
    requested selection that this table refuses to run.

    The caller decides what an unsupported cell in the selection means, and
    the sweep's answer is: it fails (see run_e2e.py). A default sweep
    (cases=None, backends=None) selects every cell of the table, so the
    unsupported ones are always PRINTED; an explicit selection that names an
    unsupported cell is a hard failure, because the caller asked for coverage
    that does not exist.
    """
    sel_cases = tuple(cases) if cases else CASES
    sel_backends = tuple(backends) if backends else BACKENDS
    for c in sel_cases:
        if c in REFERENCE_ONLY:
            # Naming a withdrawn case is a hard failure, not a quiet drop: the
            # caller asked for coverage this table does not provide, and
            # returning an empty-but-green selection would answer "this is
            # fine" to the question "did you check it" (rule 2).
            raise ValueError(
                '%s is REFERENCE_ONLY -- it is kept as reference but not gated,'
                ' so the sweep cannot run it.\n  %s' % (c, REFERENCE_ONLY[c]))
        if c not in CASES:
            raise ValueError('unknown case %r -- known cases: %s' % (c, ', '.join(CASES)))
    for b in sel_backends:
        if b not in BACKENDS:
            raise ValueError('unknown backend %r -- known backends: %s'
                             % (b, ', '.join(BACKENDS)))
    runnable, unsupported = [], []
    for c in sel_cases:
        for b in sel_backends:
            if is_supported(c, b):
                runnable.append((c, b))
            else:
                unsupported.append((c, b, UNSUPPORTED[(c, b)]))
    return runnable, unsupported


def coverage_report(runnable, declared_unsupported, selection_label):
    """The lines a run prints BEFORE it starts, so its own output states
    exactly what it is about to cover. Requirement zero of this sweep: a green
    result must never be readable as broader than it is."""
    total = len(CASES) * len(BACKENDS)
    selected = len(runnable) + len(declared_unsupported)
    lines = [
        'selection: %s' % selection_label,
        'coverage : %d of %d cells in the full %d case x %d backend table'
        % (selected, total, len(CASES), len(BACKENDS)),
        'cases    (%d): %s' % (len(CASES), ', '.join(CASES)),
        'backends (%d): %s' % (len(BACKENDS), ', '.join(BACKENDS)),
        'will RUN %d cell(s):' % len(runnable),
    ]
    for c, b in runnable:
        lines.append('  %-16s %-13s artifacts=%-7s gate: %s'
                     % (c, b, '+'.join(ARTIFACTS[b]), gate_description(c)))
    lines.append('declared UNSUPPORTED, %d cell(s) -- declared, not absent:'
                 % len(declared_unsupported))
    for c, b, reason in declared_unsupported:
        lines.append('  %-16s %-13s %s' % (c, b, reason))
    chosen = set(runnable) | set((c, b) for c, b, _ in declared_unsupported)
    not_selected = [(c, b) for c in CASES for b in BACKENDS if (c, b) not in chosen]
    lines.append('NOT in this selection, %d cell(s)%s'
                 % (len(not_selected),
                    (': ' + ', '.join('%s x %s' % cb for cb in not_selected))
                    if not_selected else ''))
    # Printed on EVERY run, selection or not: a withdrawn case is coverage this
    # sweep does not have, and the output has to say so out loud or a green
    # result quietly widens over time.
    lines.append('REFERENCE ONLY, withdrawn from the sweep, %d case(s) x %d '
                 'backend(s) = %d cell(s) NOT covered by any result below:'
                 % (len(REFERENCE_ONLY), len(BACKENDS),
                    len(REFERENCE_ONLY) * len(BACKENDS)))
    for c, reason in sorted(REFERENCE_ONLY.items()):
        lines.append('  %-16s %s' % (c, reason))
    return lines


# --- consistency, enforced at import time -------------------------------------
# A case added to testNameList.py with no entry here would otherwise run ungated
# or crash mid-sweep. Fail at import, naming the case.
_unknown_withheld = [c for c in REFERENCE_ONLY if c not in ALL_CASES]
if _unknown_withheld:
    raise RuntimeError('REFERENCE_ONLY names cases testNameList.py does not '
                       'have: %r' % _unknown_withheld)
# A withdrawn case must not ALSO carry a gate -- one state per case, or the
# table says both "not covered" and "covered at this bound".
_double_declared = [c for c in REFERENCE_ONLY if c in CASE_BOUND or c in GATE]
if _double_declared:
    raise RuntimeError('%r are REFERENCE_ONLY but still carry a bound or a '
                       'gate -- a case is gated or it is not' % _double_declared)
# "Kept as reference" is enforced, not assumed: the whole point of withdrawing
# rather than deleting is that the committed reference survives. If it is ever
# removed as unused, that is a hard failure naming the file.
for _c in REFERENCE_ONLY:
    _ref = os.path.join(REPO_ROOT, 'test.reference.results', _c,
                        'frt.canonical.txt')
    if not os.path.isfile(_ref):
        raise RuntimeError(
            '%s is REFERENCE_ONLY -- withdrawn from the sweep precisely so its '
            'reference is KEPT -- but %s is missing. Restore it from git; a '
            'withdrawn case with no reference is just a deleted case.'
            % (_c, _ref))

_missing_bound = [c for c in CASES if c not in CASE_BOUND]
_missing_gate = [c for c in CASES if c not in GATE]
_extra = [c for c in set(CASE_BOUND) | set(GATE) if c not in CASES]
if _missing_bound or _missing_gate or _extra:
    raise RuntimeError(
        'testsys/matrix.py is out of step with testNameList.py: '
        'missing CASE_BOUND for %r, missing GATE for %r, entries for unknown '
        'cases %r. Every gated case declares its bound and its gate here.'
        % (_missing_bound, _missing_gate, _extra))
for _c, _g in GATE.items():
    if _g not in ('abs-max', 'flip-budget'):
        raise RuntimeError('%s has unknown gate %r' % (_c, _g))
    if _g == 'abs-max' and CASE_BOUND[_c] is None:
        raise RuntimeError('%s is gated abs-max but has no bound' % _c)
    if _g == 'flip-budget' and CASE_BOUND[_c] is not None:
        raise RuntimeError('%s is gated flip-budget but also carries a scalar '
                           'bound -- one gate per case' % _c)
# UNSUPPORTED and CI_CELLS are statements about the SWEPT table, so they are
# checked against CASES. MEASURED_PEAK_RSS_GB is a record of measurements taken,
# including on a case since withdrawn, so it is checked against ALL_CASES --
# deleting a measurement because the case left the sweep would discard evidence.
for (_c, _b) in list(UNSUPPORTED) + list(CI_CELLS):
    if _c not in CASES or _b not in BACKENDS:
        raise RuntimeError('table entry for unknown or ungated cell %r x %r'
                           % (_c, _b))
for (_c, _b) in list(MEASURED_PEAK_RSS_GB):
    if _c not in ALL_CASES or _b not in BACKENDS:
        raise RuntimeError('measurement for unknown cell %r x %r' % (_c, _b))

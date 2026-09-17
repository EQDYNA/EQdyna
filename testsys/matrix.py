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

# The case list is testNameList.py's, not a second copy of it (rule 1). A case
# added there with no entry below is a hard import-time failure, not a silently
# ungated case -- see the consistency block at the bottom of this module.
CASES = tuple(_NAME_LIST)
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
    # TIGHTENED 2026-09-17 (item 37 -- same defect item 19a caught and fixed
    # for tpv36's bound, found here by inspection of this file's own
    # commented history a few lines below, not by a fresh claim): the
    # THRESHOLD these three carried was justified by "no python run has ever
    # been measured for these", which stopped being true once haruto's
    # CI_CELLS widening pass measured all three (`/usr/bin/time -v`,
    # full-length, one at a time, this dev box, 2026-09-16) AND confirmed
    # each one passes its cell via testsys.compare.compare_cell on that same
    # run. Observed: meng2023a numpy/jax 4.628301e-09/3.229082e-09,
    # meng2023cb numpy/jax 5.054473e-09/4.395842e-09, tpv29 numpy/jax
    # 9.876230e-15/1.276548e-13. Bounds set the same way tpv36's was: the
    # worst observation times test.tpv8's own headroom ratio
    # (1e-8 / 8.23e-11 = ~121.5x), then rounded to the nearest bound already
    # in use at that order of magnitude rather than an arbitrary new number.
    'test.meng2023a': 1e-6,    # worst 4.63e-09 * ~121.5 = 5.6e-7; rounds to
                               # the same 1e-6 already used for tpv10/tpv36.
    'test.meng2023cb': 1e-6,   # worst 5.05e-09 * ~121.5 = 6.1e-7; same bound.
    'test.tpv29': 1e-10,       # worst 1.28e-13 * ~121.5 = 1.6e-11; rounded UP
                               # to 1e-10 (headroom ~780x, not the smaller
                               # ~121x -- these are near-machine-epsilon
                               # diffs and 1e-11 would leave no margin for
                               # ordinary cross-run noise at that scale).
    'test.drv.a6': None,      # chaotically bistable -- flip-budget gate, below.
    # C_degen>3 (wedge-degeneration), dip=15 (tpv36's par.dip), dx=500,
    # par.term=6 -- a COARSE regression gate (rule 17 step 4), NOT a claim
    # about physical accuracy at SCEC-standings resolution (the published
    # tpv36/37 standings rank a run this coarse last and improve
    # monotonically with resolution). Measured AFTER item 19(a)'s jax-only
    # regression fix (assembleGlobalMass.py's eleshp einsum dispatch), on a
    # reference frozen fresh for this commit -- the merge's original
    # measurement predates that fix and is not reused: observed python-jax
    # 5.867293e-09, python-numpy 7.894064e-09 (both near-zero peak-slip-rate
    # components on marginally-ruptured nodes -- roundoff, not physics).
    # 1e-6 is ~127x the worst observation, the same headroom ratio
    # test.tpv8 carries (1e-8 / 8.23e-11).
    'test.tpv36': 1e-6,
    # test.tpv37: identical wedge-degeneration machinery as tpv36 (same
    # dip=15, dx=500, par.term=6, C_degen=par.dip) -- only par.tpv=37 and one
    # initial-stress-taper coefficient differ from tpv36's compset (confirmed
    # by diff: 0.001875e6 vs tpv36's 0.0005e6 in on_fault_vars[iz,ix,4]).
    # Same coarse-regression-gate reasoning as tpv36 (rule 17 step 4, NOT a
    # SCEC-standings accuracy claim). Fortran cell bit-exact against a
    # freshly-frozen reference (max|diff|=0.0). Measured fresh, this commit:
    # python-jax 5.256410e-09, python-numpy 6.929140e-09 (both near-zero
    # peak-slip-rate-type components on marginally-ruptured nodes -- roundoff,
    # not physics, same shape as tpv36's own observation). 1e-6 is ~144x the
    # worst observation, in the same headroom family as tpv8/tpv36.
    'test.tpv37': 1e-6,
}

GATE = {
    'test.tpv8': 'abs-max',
    'test.tpv104': 'abs-max',
    'test.tpv1053d': 'abs-max',
    'test.tpv10': 'abs-max',
    'test.meng2023a': 'abs-max',
    'test.meng2023cb': 'abs-max',
    'test.tpv29': 'abs-max',
    'test.tpv36': 'abs-max',
    'test.tpv37': 'abs-max',
    # test.drv.a6 (C_elastic==0 viscoplastic, friclaw==4, fractal-rough, long
    # duration) has genuinely bistable rupture arrivals, so a scalar max-abs
    # bound cannot distinguish "a few hundred marginal nodes flipped" from
    # "everything drifted a little". Gated on bulk agreement PLUS an explicit
    # flip budget instead.
    'test.drv.a6': 'flip-budget',
}

# test.drv.a6's flip-budget gate. These numbers are measured, not chosen;
# testsys/parity/evidence_drv_a6_chaos.py regenerates them and they must not be
# changed without re-running it fresh (rule 4).
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

# Cells CI runs, and the measured reason the rest are left out. Historically
# (one job) a GitHub ubuntu-22.04 runner's 7 GB was shared by the Fortran
# build and every python/jax cell at once; exceeding it produced a SIGTERM
# with no output (exit 143), which is a resource kill and not a test result.
# As of 2026-09-16 the workflow is matrixed into parallel jobs, each with its
# own 7 GB runner -- see the CI_CELLS comment below for what that changed and
# what it did not. This is a declared, measured decision with the numbers
# next to it -- not an env var whose value nobody can audit.
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
    # Measured 2026-09-16 (haruto, CI-widening pass), `/usr/bin/time -v`,
    # full-length, one at a time, against the CURRENT tree (a first attempt
    # ran against a tree that still had mira's since-reverted C_degen merge
    # and was discarded and re-measured fresh -- rule 4, a stale-tree run is
    # not evidence). Every one of these also PASSED its bound on this same
    # run (testsys.compare.compare_cell), not just produced a number:
    # meng2023a numpy/jax 4.628301e-09/3.229082e-09, meng2023cb numpy/jax
    # 5.054473e-09/4.395842e-09, tpv29 numpy/jax 9.876230e-15/1.276548e-13,
    # all against THRESHOLD=1e-3.
    ('test.meng2023a', 'python-numpy'): 2.545,    # 313.4 s wall
    ('test.meng2023a', 'python-jax'): 4.025,      # 70.7 s wall
    ('test.meng2023cb', 'python-numpy'): 2.543,   # 419.7 s wall
    ('test.meng2023cb', 'python-jax'): 4.020,     # 65.1 s wall
    ('test.tpv29', 'python-numpy'): 2.705,   # 1711.0 s wall -- by far the
                                              # most expensive cell in the
                                              # table; isolated in its own CI
                                              # job (e2e-ci-python-tpv29) so
                                              # nothing else queues behind it
                                              # and it never queues behind
                                              # anything else.
    ('test.tpv29', 'python-jax'): 4.052,     # 198.0 s wall
}
CI_RUNNER_RAM_GB = 7.0
# WIDENED 2026-09-16 (haruto, per owner instruction). The reasoning that
# excluded tpv10/tpv104/tpv1053d x python-jax was STALE: it was written when
# ONE job held the Fortran build, mpirun, AND every python/jax cell at once,
# so their 3.2-4.0 GB looked like it was competing with everything else for
# the same 7 GB. Once the workflow is matrixed (this repo, 2026-09-16) each
# job gets its own 7 GB, AND run_e2e.py's own core-budget allocator already
# forces every cell in a job to run ONE AT A TIME on a small runner (verified
# on two real CI runs: job logs show cells finishing back-to-back, wall clock
# == sum of per-cell seconds, never a "running N cells concurrently" line) --
# so a job's peak memory is the single LARGEST cell it holds, not a sum, and
# always was, even in the original one-job workflow. On the measured numbers
# alone only test.drv.a6 x python-jax (9.57 GB) actually exceeds a 7 GB
# runner; the rest below were excluded by a constraint the matrix removed.
CI_CELLS = (
    tuple((c, 'fortran') for c in CASES)
    + (('test.tpv8', 'python-numpy'), ('test.tpv8', 'python-jax'),
       ('test.tpv10', 'python-jax'), ('test.tpv104', 'python-jax'),
       ('test.tpv1053d', 'python-jax'),
       ('test.meng2023a', 'python-numpy'), ('test.meng2023a', 'python-jax'),
       ('test.meng2023cb', 'python-numpy'), ('test.meng2023cb', 'python-jax'),
       ('test.tpv29', 'python-numpy'), ('test.tpv29', 'python-jax'))
)
# WHAT THIS LIST LEAVES OUT, AND WHY -- said here rather than implied. Every
# exclusion below is a MEASURED-or-genuinely-unmeasured decision; none is a
# case that fails.
#   * test.drv.a6 x python-jax: 9.57 GB measured -- the one real exclusion.
#   * test.drv.a6 x python-numpy: never measured. Not chased this round (owner
#     instruction, 2026-09-16): drv.a6 is the largest/most expensive case in
#     the table by a wide margin and its jax column already leaves a 7 GB
#     runner with zero margin; a larger runner is a separate decision.
#   * test.tpv10/tpv104/tpv1053d x python-NUMPY (their jax columns are now IN
#     CI_CELLS, see above): never measured. The rule stays literal -- a cell
#     without a measured RSS does not go into CI, and "numpy is probably in
#     the same ballpark as jax" is exactly the kind of guess rule 6 forbids.
#     Cheapest next step for whoever widens further: measure these three the
#     same way as the entries above.
#   * test.tpv36 x python-numpy/python-jax: same reason -- peak RSS not yet
#     measured on CI's runner. tpv36's fortran cell joined CI_CELLS
#     automatically (every case's fortran cell does, unconditionally, by this
#     tuple's own construction above) the moment the case was re-gated
#     (item 19a, after the jax-only regression that had reverted it was
#     found and fixed -- see pathway_forward.md).
#   * test.tpv37 x python-numpy/python-jax: same reason as tpv36 -- peak RSS
#     not yet measured on CI's runner. tpv37's fortran cell joined CI_CELLS
#     the same automatic way the moment it was gated (item 19a, second half).
# A green CI run therefore means 21 of 30 cells, and says so.


def is_supported(case, backend):
    return (case, backend) not in UNSUPPORTED



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
    return lines


# --- consistency, enforced at import time -------------------------------------
# A case added to testNameList.py with no entry here would otherwise run ungated
# or crash mid-sweep. Fail at import, naming the case.
# Every gated case must have a committed reference. Absent = broken checkout,
# and a checkout that cannot compare must not be able to report a pass.
for _c in CASES:
    _ref = os.path.join(REPO_ROOT, 'test.reference.results', _c,
                        'frt.canonical.txt')
    if not os.path.isfile(_ref):
        raise RuntimeError(
            '%s is gated but %s is missing -- restore it from git. A gated case '
            'with no reference cannot be compared, and "could not compare" must '
            'never read as "passed" (rule 2).' % (_c, _ref))

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
for (_c, _b) in list(UNSUPPORTED) + list(CI_CELLS):
    if _c not in CASES or _b not in BACKENDS:
        raise RuntimeError('table entry for unknown or ungated cell %r x %r'
                           % (_c, _b))
for (_c, _b) in list(MEASURED_PEAK_RSS_GB):
    if _c not in CASES or _b not in BACKENDS:
        raise RuntimeError('measurement for unknown cell %r x %r' % (_c, _b))

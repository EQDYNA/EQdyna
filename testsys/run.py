#! /usr/bin/env python3
"""
Single entry point for EQdyna's tiered test system (PROJECT_RULES.md rule 3).

    python3 testsys/run.py unit          # fast pure-python unit tests (pytest, no MPI/Fortran)
    python3 testsys/run.py regression    # one guard per past incident (rule 10)
    python3 testsys/run.py e2e           # THE test: the everyday case x backend sweep vs test.reference.results/, at matrix.GATE_TERM_S (rule 7)
    python3 testsys/run.py e2e-ci        # the same sweep, restricted to matrix.CI_CELLS (2026-09-23: one case x 3 backends, portability smoke)
    python3 testsys/run.py release       # THE test, over every supported cell (everyday cells + matrix.RELEASE_ONLY): the release gate, writes docs/evidence/sweep-<sha>/summary.json
    python3 testsys/run.py e2e-full      # SCEC cases at spec dx/term, 16 ranks, report-only (opt-in; hours)
    python3 testsys/run.py gpu           # the sweep's python-jax column on CUDA (one cell; needs a GPU)
    python3 testsys/run.py perf          # pinned single-core Fortran/NumPy/JAX timing, ratio-guarded
    python3 testsys/run.py all           # unit + regression + e2e, in that order (default; perf is opt-in, not in "all" -- it needs a Fortran build a fresh checkout does not have yet)

There is ONE test here -- e2e -- and backend is an axis of it, not a tier.
There is also ONE term (2026-09-23 owner decision, superseding a same-day
earlier two-term design): every cell of every selection -- `e2e`, `e2e-ci`,
`release` alike -- runs at matrix.GATE_TERM_S (5 s), applied through the ONE
existing override path regardless of a case's own committed par.term.
`release` differs from `e2e` only in which CELLS run: every supported cell
(everyday cells + matrix.RELEASE_ONLY, run_e2e.py's `--release` selector), not
a different term -- and it is the tier that gates a release (PROJECT_RULES
rule 15/16). CI runs e2e-ci ONLY; `release` is a human- or conductor-scheduled
tier, the same as the old `run.py all`'s e2e slice was.

Two tiers that used to sit beside it are gone, for the same reason each time:
they were a SECOND way of asking a question e2e already answers, and a second
way of asking costs a second implementation to keep honest.
  * `accept` -- the same solver against the same references with its own
    comparison implementation over a shorter case list.
  * `parity` -- Python vs a Fortran serial oracle, per step, via the
    eqdyna-pydump build. It answered "at WHICH STEP does it diverge?", which
    is a debugging question, not a gating one; it root-caused the missing
    rsfNucleation branch and that work is done. It went stale the moment the
    Python tree was restructured -- its two test files were named for a
    `standalone/` package that no longer exists -- and keeping a stale tier
    green costs more than re-deriving a step dump the next time one is needed.
    Removed with it: make_fixtures.py, the pydump fixtures, and src/pydump*
    (the Fortran link-time seam they existed to drive).

Prints a per-test SUCCESS/FAIL line (from pytest or from each regression
script's own banner), a per-tier SUMMARY line, and exits non-zero if
anything in the requested scope failed. "The script ran" and "the script
passed" are always two different questions here (rule 3).
"""
import glob
import os
import subprocess
import sys

TESTSYS = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(TESTSYS)
TIERS = ('unit', 'regression', 'e2e')

# Line-buffered: redirected to a log (CI pipes this through `tee`), Python
# block-buffers this script's own prints while the children write straight to
# the fd, so the tier banners and the SUMMARY would land out of order, and
# would vanish entirely if the job were killed.
sys.stdout.reconfigure(line_buffering=True)


def run_unit():
    d = os.path.join(TESTSYS, 'unit')
    print('\n==== testsys: unit ====')
    return subprocess.call([sys.executable, '-m', 'pytest', '-v', d], cwd=REPO_ROOT)


def run_regression():
    """Each regression/test_*.py is a standalone script with its own
    SUCCESS/FAIL banner and sys.exit (matching the shape of the incident
    it guards) rather than pytest functions -- run each one and gate on
    its exit code."""
    print('\n==== testsys: regression ====')
    scripts = sorted(glob.glob(os.path.join(TESTSYS, 'regression', 'test_*.py')))
    if not scripts:
        print('regression: FAIL - no regression scripts found (misconfigured testsys/)')
        return 1
    overall = 0
    for script in scripts:
        name = os.path.basename(script)
        print(f'-- {name} --')
        rc = subprocess.call([sys.executable, script], cwd=REPO_ROOT)
        print(f'{"SUCCESS" if rc == 0 else "FAIL"} {name} (exit {rc})')
        overall = overall or rc
    return overall


def _e2e(*args):
    return subprocess.call([sys.executable, os.path.join(TESTSYS, 'e2e', 'run_e2e.py')]
                           + list(args), cwd=REPO_ROOT)


def run_e2e():
    """THE sweep, EVERYDAY selection (run_e2e.py's default, no flags): every
    case minus matrix.RELEASE_ONLY, at matrix.GATE_TERM_S regardless of its
    own committed par.term."""
    print('\n==== testsys: e2e ====')
    return _e2e()


def run_e2e_ci():
    """The sweep restricted to matrix.CI_CELLS -- 2026-09-23 (owner-approved
    test-methodology change): one case (test.tpv8), all three backend
    implementations, at the gate term. Its job is portability (clean
    checkout, fresh deps, a different MPI from this box's), not physics
    coverage -- that job belongs to `e2e` and `release` now. This is a
    declared selection, not a magic env var: the run prints which cells it
    covered and which it did not."""
    print('\n==== testsys: e2e-ci ====')
    return _e2e('--ci')


def run_release():
    """THE sweep, RELEASE selection (run_e2e.py --release): every supported
    cell (everyday cells + matrix.RELEASE_ONLY), at the SAME matrix.GATE_TERM_S
    as every other selection -- this tier widens which CELLS run, not the
    term. This is the release gate (PROJECT_RULES rule 15/16) -- the tier
    `run.py all`'s e2e slice used to be before CI stopped running the wider
    sweep. Run over the full default selection (no --cases/--backends), it
    also writes docs/evidence/sweep-<shortsha>/summary.json
    (run_e2e.write_release_evidence).

    `run.py unit regression release` in ONE invocation is SAFE: checked
    2026-09-23 (rule-24 tree_clean fix) that no regression script under
    testsys/regression/ writes to a real tracked repo path -- every script
    that touches docs/perf_ledger.jsonl-shaped data does so inside its own
    tempfile.mkdtemp() sandbox or reads the committed ledger read-only
    (test_perf_ledger.py); test_sweep_speed_2026_09_23.py imports run_e2e
    but calls schedule_order() directly, never main(), so it runs no sweep.
    run_e2e.py's tree_clean is captured at the very start of ITS OWN
    process (run_e2e.capture_start_tree_state, called first thing in
    main()), so even if a future regression script broke this invariant by
    writing to a real tracked path, `release`'s reading would still be
    correct -- it would just, correctly, read as dirty. This tier is not
    required to run standalone; nothing before it in the same invocation
    is known to dirty the tree, and if that ever changes the release
    evidence will say so rather than silently passing."""
    print('\n==== testsys: release ====')
    return _e2e('--release')


def run_gpu():
    """The sweep's python-jax column on CUDA. A selection of the one sweep,
    not a tier of its own: JAX_PLATFORMS=cuda, same cases, same references,
    same per-case bounds. Fails loudly if there is no CUDA device -- "no GPU
    here" and "the GPU column passed" must not share an exit code (rule 2)."""
    print('\n==== testsys: gpu ====')
    try:
        import jax
    except Exception as exc:                        # noqa: BLE001
        print('FAIL gpu: jax is not importable (%s) -- this selection cannot '
              'be gated on this machine' % exc)
        return 1
    devs = [d for d in jax.devices() if 'cuda' in str(d).lower() or 'gpu' in str(d).lower()]
    if not devs:
        print('FAIL gpu: no CUDA device/plugin visible to jax (pip install '
              '"jax[cuda12]" on a GPU box). Reporting this as a failure, not a '
              'skip: a green line here would claim GPU coverage that did not '
              'happen.')
        return 1
    return _e2e('--backends', 'python-jax', '--cases', 'test.tpv8',
                '--device', 'cuda')


def run_e2e_full():
    print('\n==== testsys: e2e-full ====')
    return subprocess.call(
        [sys.executable, os.path.join(TESTSYS, 'e2e', 'run_e2e.py'), '--full'],
        cwd=REPO_ROOT)


def run_scaling():
    print('\n==== testsys: scaling ====')
    return subprocess.call([sys.executable, os.path.join(TESTSYS, 'perf', 'run_scaling.py')], cwd=REPO_ROOT)


def run_perf():
    print('\n==== testsys: perf ====')
    return subprocess.call([sys.executable, os.path.join(TESTSYS, 'perf', 'run_perf.py')],
                            cwd=REPO_ROOT)


RUNNERS = {'unit': run_unit, 'regression': run_regression, 'e2e': run_e2e,
           'e2e-ci': run_e2e_ci, 'release': run_release, 'perf': run_perf,
           'gpu': run_gpu, 'scaling': run_scaling, 'e2e-full': run_e2e_full}
# 'all' stays unit+regression+e2e only (TIERS below, `e2e` the everyday
# selection) -- perf requires a Fortran build and a generated baseline that a
# fresh checkout does not have; it is opt-in, invoked by name, not swept into
# 'all'. e2e-ci and gpu are SELECTIONS of e2e, not tiers of their own;
# `release` is the wider (everyday + RELEASE_ONLY) selection of the same one
# sweep, at the SAME GATE_TERM_S, likewise not a second implementation. 'all'
# runs the everyday sweep; CI runs e2e-ci (2026-09-23: a one-case portability
# smoke, not physics coverage); `release` is the human-/conductor-scheduled
# release gate, whose narrower or wider coverage the sweep itself prints
# either way.
# e2e-full additionally needs EQDYNA_FULL_LAUNCH=yes-hours (see
# testsys/e2e/run_e2e_full.py) -- spec-resolution SCEC runs are hours long
# and user-scheduled, never automatic.
OPTIONAL_TIERS = ('e2e-ci', 'release', 'perf', 'gpu', 'scaling', 'e2e-full')


def main(argv):
    # Several tiers may be named in one invocation, in the order given --
    # that is how CI asks for "unit regression e2e-ci" without a second
    # entry point and without an env var deciding coverage behind its back.
    requested = argv[1:] or ['all']
    unknown = [t for t in requested if t not in TIERS + OPTIONAL_TIERS + ('all',)]
    if unknown:
        print(f'unknown tier(s) {unknown}')
        print(f'usage: python3 testsys/run.py [{"|".join(TIERS + OPTIONAL_TIERS)}|all] ...')
        return 2

    selected = []
    for tier in requested:
        for t in (TIERS if tier == 'all' else (tier,)):
            if t not in selected:
                selected.append(t)

    results = {tier: RUNNERS[tier]() for tier in selected}

    print('\n==== testsys: SUMMARY ====')
    print('tiers run: %s' % ', '.join(selected))
    overall = 0
    for tier in selected:
        rc = results[tier]
        print(f'{"SUCCESS" if rc == 0 else "FAIL"} {tier} (exit {rc})')
        overall = overall or rc

    return 1 if overall else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))

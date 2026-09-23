#! /usr/bin/env python3
"""
Regression guard: the 2026-09-23 sweep-speed changes (rules 2, 10).

RELEASE_ONLY (matrix.py, item 1). test.tpv36 x python-numpy and test.tpv37 x
python-numpy are SUPPORTED and PASSING, held out of the EVERYDAY sweep
(run.py e2e / run_e2e.py's default selection at --term gate) for wall-clock
cost only, and restored at --term full (run.py release, rule 24). The defect
shape this guards against: a third cell added to RELEASE_ONLY with no
matching update here goes RED (a silent widening of what "everyday" no
longer covers must not pass quietly), and a RELEASE_ONLY cell that is also
declared UNSUPPORTED goes RED (release-only means "supported, cost-
deferred", never "does not work").

Cheap (rule 9): imports + plain-data assertions. No solver, no MPI, no I/O
beyond reading matrix.py. Exits non-zero on any failure.

NOTE: this file grows a second section (item 2, longest-first scheduling) in
a follow-up change to the same session -- see its own commit for that half.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from testsys import matrix  # noqa: E402


# --------------------------------------------------------------------------
# item 1: RELEASE_ONLY
# --------------------------------------------------------------------------
EXPECTED_RELEASE_ONLY = {
    ('test.tpv36', 'python-numpy'),
    ('test.tpv37', 'python-numpy'),
}


def check_release_only_is_exactly_the_two_flagged_cells():
    """Pins the CONTENT of the flag: a third cell added here with no update
    to EXPECTED_RELEASE_ONLY above goes RED, and so does one of the two
    being silently dropped."""
    actual = set(matrix.RELEASE_ONLY)
    assert actual == EXPECTED_RELEASE_ONLY, (
        'matrix.RELEASE_ONLY is %r, expected exactly %r -- if a cell was '
        'deliberately added or removed from the release-only flag, update '
        'EXPECTED_RELEASE_ONLY in this guard as part of that same change'
        % (actual, EXPECTED_RELEASE_ONLY))
    print('  PASS  matrix.RELEASE_ONLY == %r' % sorted(EXPECTED_RELEASE_ONLY))


def check_release_only_cells_are_never_declared_unsupported():
    """release-only means "supported, cost-deferred to the release tier" --
    never "does not work". matrix.py's own import-time consistency check
    already refuses this combination; this restates the contract at the
    level a reader of this guard will look for it."""
    overlap = set(matrix.RELEASE_ONLY) & set(matrix.UNSUPPORTED)
    assert overlap == set(), (
        '%r is both RELEASE_ONLY and UNSUPPORTED -- a cell cannot honestly '
        'claim both' % overlap)
    print('  PASS  no RELEASE_ONLY cell is declared UNSUPPORTED')


def check_release_selection_equals_everyday_selection_plus_release_only():
    """The two selections run_e2e.py's default (non-explicit) path chooses
    between: --term gate (everyday) must be exactly --term full (release)
    minus RELEASE_ONLY, and vice versa. This is the mutation-tested
    boundary: move a cell across it by hand (below) and confirm it flips."""
    release_runnable, release_unsupported = matrix.cells()
    everyday_runnable, everyday_unsupported, everyday_release_only = \
        matrix.everyday_cells()

    assert set(release_runnable) - set(everyday_runnable) == EXPECTED_RELEASE_ONLY, (
        'release selection minus everyday selection must be exactly the '
        'flagged cells')
    assert set(everyday_runnable) & EXPECTED_RELEASE_ONLY == set(), (
        'everyday selection must exclude every RELEASE_ONLY cell')
    assert set((c, b) for c, b, _ in everyday_release_only) == EXPECTED_RELEASE_ONLY
    # unsupported cells are identical on both sides -- RELEASE_ONLY only
    # ever moves a cell between "runnable" and "held back", never touches
    # the unsupported set.
    assert release_unsupported == everyday_unsupported
    print('  PASS  release == everyday + %d release-only cell(s), '
          'unsupported set unchanged' % len(EXPECTED_RELEASE_ONLY))


def check_mutation_moving_a_cell_across_the_boundary_flips_the_guard():
    """Take one real, currently-everyday cell, move it across the boundary
    by hand (a synthetic RELEASE_ONLY dict, not matrix.py's real one), and
    confirm the equality above flips red, then confirm the untouched real
    table still reports the original (unflipped) relationship -- proving
    this guard actually distinguishes the two states rather than being
    trivially true either way."""
    real_release_only = dict(matrix.RELEASE_ONLY)
    probe_cell = ('test.tpv8', 'python-jax')  # a real, everyday, non-flagged cell
    assert probe_cell not in real_release_only, (
        '%r is already RELEASE_ONLY -- pick a different probe cell for this '
        'mutation check' % (probe_cell,))

    mutated = dict(real_release_only)
    mutated[probe_cell] = 'synthetic, test-only'
    try:
        matrix.RELEASE_ONLY = mutated
        _, _, release_only_after = matrix.everyday_cells()
        moved = set((c, b) for c, b, _ in release_only_after)
        assert probe_cell in moved, (
            'moving %r into RELEASE_ONLY did not remove it from the '
            'everyday selection -- the boundary check above would not '
            'catch a real regression here' % (probe_cell,))
        assert moved != EXPECTED_RELEASE_ONLY, (
            'the mutated 3-cell release-only set was not distinguishable '
            'from the real 2-cell one')
    finally:
        matrix.RELEASE_ONLY = real_release_only

    # restored: the real table reports exactly the original relationship again
    _, _, release_only_restored = matrix.everyday_cells()
    assert set((c, b) for c, b, _ in release_only_restored) == EXPECTED_RELEASE_ONLY
    print('  PASS  moving %r across the boundary flips the selection and '
          'restoring the table restores it (mutation-tested both ways)'
          % (probe_cell,))


def main():
    print('Regression guard: 2026-09-23 sweep-speed changes (RELEASE_ONLY)')
    checks = [
        check_release_only_is_exactly_the_two_flagged_cells,
        check_release_only_cells_are_never_declared_unsupported,
        check_release_selection_equals_everyday_selection_plus_release_only,
        check_mutation_moving_a_cell_across_the_boundary_flips_the_guard,
    ]
    failures = []
    for c in checks:
        try:
            c()
        except AssertionError as e:
            failures.append('%s: %s' % (c.__name__, e))
            print('  FAIL  %s: %s' % (c.__name__, e))
    if failures:
        print('\nFAIL test_sweep_speed_2026_09_23 (%d check(s))' % len(failures))
        return 1
    print('\nSUCCESS test_sweep_speed_2026_09_23')
    return 0


if __name__ == '__main__':
    sys.exit(main())

#! /usr/bin/env python3
"""
Regression guard: the mechanical profile guard (testsys/profile_schema.py's
`validate()`/`validate_run_dir()`) must go RED on the three ways a per-rank
profile.rank<r>.json can lie, and stay GREEN on real emitter output from
every backend that writes one (fortran, python-numpy, python-jax,
python-jax-mpi).

WHY THIS MATTERS (docs/run_profile.md, CLAUDE.md "things that will bite
you"): `writeCompTime` initialised to 0 and read from nowhere let a 511x
error in `compTimeInSeconds(2)` survive completely undetected, because
nothing ever summed that bucket against anything. This guard is the thing
that now sums it, on every run, for every backend -- and this file is the
proof that the guard itself actually catches the shape of error it exists
for, not merely that it runs without crashing on tidy input.

FIXTURES (testsys/regression/fixtures/profile_guard/<backend>/): REAL
`profile.rank<r>.json` output, copied verbatim (unedited) from one real
`testsys/e2e/run_e2e.py --cases test.tpv8 --backends <all four>` sweep run
on this box, 2026-09-23 (this session; see docs/notes/NOTES_profile_guard.md). Not
synthetic: this is what the emitters actually wrote, four fortran ranks
(4-rank test.tpv8), one python-numpy rank, one python-jax rank, four
python-jax-mpi ranks. The GREEN checks below assert the guard accepts these
UNCHANGED; the RED checks mutate a throwaway in-memory copy or a throwaway
tempdir copy and never touch the committed fixture files.

Cheap (rule 9): no solver launch, no subprocess -- pure file/dict validation
against already-committed fixtures. Milliseconds. Exits non-zero on any
failure.
"""
import copy
import json
import os
import shutil
import sys
import tempfile

TESTSYS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO_ROOT = os.path.dirname(TESTSYS)
sys.path.insert(0, TESTSYS)
import profile_schema  # noqa: E402

FIXTURES = os.path.join(TESTSYS, 'regression', 'fixtures', 'profile_guard')

# (backend, dirname, expected rank count) -- the real ranks each fixture
# subdirectory was captured at (see this file's own docstring).
BACKEND_FIXTURES = (
    ('fortran', 'fortran', 4),
    ('python-numpy', 'python-numpy', 1),
    ('python-jax', 'python-jax', 1),
    ('python-jax-mpi', 'python-jax-mpi', 4),
)


def _load_one_row(backend_dir):
    """The rank0 row of a fixture, as a fresh dict (safe to mutate)."""
    path = os.path.join(FIXTURES, backend_dir, 'profile.rank0.json')
    with open(path) as f:
        return json.load(f)


# --------------------------------------------------------------------------
# GREEN: real emitter output from all four backends must validate as-is
# --------------------------------------------------------------------------
def _check_green_one_backend(backend, backend_dir, expected_ranks):
    d = os.path.join(FIXTURES, backend_dir)
    rows = profile_schema.validate_run_dir(d)
    assert len(rows) == expected_ranks, (
        '%s fixture at %r: validate_run_dir returned %d row(s), expected '
        '%d (the real rank count this fixture was captured at)'
        % (backend, d, len(rows), expected_ranks))
    for row in rows:
        assert row['backend'] == backend, (
            '%s fixture row claims backend=%r' % (backend, row['backend']))
        assert row['schema'] == profile_schema.SCHEMA_ID


def check_green_on_real_fortran_output():
    _check_green_one_backend(*BACKEND_FIXTURES[0])


def check_green_on_real_numpy_output():
    _check_green_one_backend(*BACKEND_FIXTURES[1])


def check_green_on_real_jax_output():
    _check_green_one_backend(*BACKEND_FIXTURES[2])


def check_green_on_real_jaxmpi_output():
    _check_green_one_backend(*BACKEND_FIXTURES[3])


# --------------------------------------------------------------------------
# RED 1: a missing rank file
# --------------------------------------------------------------------------
def check_red_on_missing_rank_file():
    """Copy the 4-rank fortran fixture to a tempdir, delete rank2's file,
    and confirm validate_run_dir names exactly rank 2 as missing -- not a
    generic failure, and not a silent pass with 3 of 4 ranks."""
    src = os.path.join(FIXTURES, 'fortran')
    with tempfile.TemporaryDirectory() as tmp:
        for name in os.listdir(src):
            shutil.copy(os.path.join(src, name), os.path.join(tmp, name))
        os.remove(os.path.join(tmp, 'profile.rank2.json'))
        try:
            profile_schema.validate_run_dir(tmp)
            raise AssertionError(
                'validate_run_dir accepted a run tree missing rank 2 of 4 '
                '-- a launch that started fewer real workers than declared '
                'must not be able to look complete')
        except ValueError as exc:
            assert 'missing rank file' in str(exc) and '[2]' in str(exc), (
                'wrong diagnosis for a missing rank file: %s' % exc)


# --------------------------------------------------------------------------
# RED 2: a malformed / missing field
# --------------------------------------------------------------------------
def check_red_on_missing_field():
    """Delete a required field from a real row (in memory) and confirm
    validate() raises naming that field -- not a KeyError leaking out of the
    guard itself, which would look like a crash rather than a diagnosis."""
    row = _load_one_row('python-numpy')
    del row['total_s']
    try:
        profile_schema.validate(row)
        raise AssertionError(
            'validate() accepted a row missing required field total_s')
    except ValueError as exc:
        assert 'total_s' in str(exc), (
            'wrong diagnosis for a missing field: %s' % exc)


def check_red_on_malformed_field_type():
    """A bucket carrying a string instead of a number -- the shape of a
    write that landed a stray unit or format artefact -- must raise, not
    coerce."""
    row = _load_one_row('python-jax')
    row['buckets_s']['element'] = 'NaN'
    try:
        profile_schema.validate(row)
        raise AssertionError(
            'validate() accepted buckets_s["element"] = "NaN" (a string, '
            'not a number)')
    except ValueError as exc:
        assert 'element' in str(exc), (
            'wrong diagnosis for a malformed bucket type: %s' % exc)


# --------------------------------------------------------------------------
# RED 3: a bucket inflated 511x -- the compTimeInSeconds(2) incident shape
# --------------------------------------------------------------------------
def check_red_on_bucket_inflated_511x():
    """Multiply one real bucket by 511 (the measured ratio of the historical
    compTimeInSeconds(2) incident: 0.0001 s before the fix, 0.0511 s after --
    see docs/run_profile.md), leaving unaccounted_s exactly as the emitter
    wrote it. The independent-remainder check must catch this on its own,
    with no help from the caller recomputing anything."""
    row = _load_one_row('fortran')
    row['buckets_s']['element'] = row['buckets_s']['element'] * 511
    try:
        profile_schema.validate(row)
        raise AssertionError(
            'validate() accepted a bucket inflated 511x while unaccounted_s '
            'was left at its original (now-stale) value -- this is exactly '
            'the compTimeInSeconds(2) incident shape')
    except ValueError as exc:
        assert 'unaccounted_s' in str(exc), (
            'wrong diagnosis for a 511x-inflated bucket: %s' % exc)


def check_red_on_bucket_inflated_511x_even_with_unaccounted_recomputed():
    """Same mutation, but unaccounted_s IS recomputed consistently with the
    inflated bucket (so check 1, the float-precision remainder check, is
    satisfied) -- this isolates check 2, the wide 5%-of-total_s tolerance,
    and confirms it alone still catches a x511 bucket."""
    row = _load_one_row('fortran')
    row['buckets_s']['element'] = row['buckets_s']['element'] * 511
    bucket_sum = sum(row['buckets_s'].values())
    row['unaccounted_s'] = row['total_s'] - bucket_sum
    try:
        profile_schema.validate(row)
        raise AssertionError(
            'validate() accepted a bucket inflated 511x even after '
            'unaccounted_s was recomputed to match it -- the wide '
            'SUM_TOLERANCE check must catch this on its own')
    except ValueError as exc:
        assert 'exceeds' in str(exc), (
            'wrong diagnosis for a x511 bucket past the wide tolerance: %s'
            % exc)


def check_sum_floor_admits_short_run_overhead_but_not_more():
    """SUM_FLOOR_S (2026-09-24): the gap may be up to max(5% of total_s,
    2 s). The master CI failure it fixes was 0.80 s of 15.40 s (5.2%) on
    test.tpv8 x python-jax. That must now pass. A 2.5 s gap on the same 15 s
    run, and a 10 s gap on a 190 s run, must still fail."""
    b = {k: 0.0 for k in profile_schema.BUCKET_KEYS}
    def verdict(total, gap):
        bb = dict(b); bb['element'] = total - gap
        try:
            profile_schema.check_buckets('floor-check', total, bb, gap)
            return True
        except ValueError:
            return False
    assert verdict(15.40, 0.80), 'the CI case (0.80 s of 15.40 s) must pass'
    assert not verdict(15.0, 2.5), 'a 2.5 s gap on a 15 s run must fail'
    assert verdict(190.0, 9.0), '9 s of 190 s (4.7%) must pass'
    assert not verdict(190.0, 10.0), '10 s of 190 s (5.3%) must fail'


def main():
    checks = [check_green_on_real_fortran_output,
              check_green_on_real_numpy_output,
              check_green_on_real_jax_output,
              check_green_on_real_jaxmpi_output,
              check_red_on_missing_rank_file,
              check_red_on_missing_field,
              check_red_on_malformed_field_type,
              check_red_on_bucket_inflated_511x,
              check_red_on_bucket_inflated_511x_even_with_unaccounted_recomputed,
              check_sum_floor_admits_short_run_overhead_but_not_more]
    failures = []
    for c in checks:
        try:
            c()
            print('  PASS  %s' % c.__name__)
        except AssertionError as e:
            failures.append('%s: %s' % (c.__name__, e))
            print('  FAIL  %s: %s' % (c.__name__, e))
    if failures:
        print('\nFAIL test_profile_guard (%d check(s))' % len(failures))
        return 1
    print('\nSUCCESS test_profile_guard (%d checks)' % len(checks))
    return 0


if __name__ == '__main__':
    sys.exit(main())

#! /usr/bin/env python3
"""Behavioural regression guard for board row 112.i (2026-09-24):
run_e2e.cell_cost/profile_ranks used `max(1, <configured ranks>)` for the
fortran and python-jax-mpi backends, which silently turns a misconfigured
non-positive matrix.FORTRAN_RANKS/PY_MPI_RANKS entry into a plausible "1
rank" cell instead of naming the bad config (rule 2: no silent fallback).

Drives the real functions with a monkeypatched matrix table -- no case
build, no solver run -- so this file is fast and machine-independent.
`matrix.JAX_MEASURED_CORES`'s own max(1, ceil(...)) is untouched by this
fix (a different defect shape: a measured-cores CEILING, not a configured
rank count that could be silently wrong) and is not exercised here.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TESTSYS = os.path.dirname(HERE)
ROOT = os.path.dirname(TESTSYS)
E2E = os.path.join(TESTSYS, "e2e")
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
sys.path.insert(0, E2E)

import run_e2e  # noqa: E402
from testsys import matrix  # noqa: E402

FAILURES = []


def check(ok, what):
    print("%s -- %s" % ("PASS" if ok else "FAIL", what))
    if not ok:
        FAILURES.append(what)


def raises(fn, exc_types, what, must_contain=()):
    try:
        fn()
    except exc_types as e:
        msg = str(e)
        absent = [f for f in must_contain if f not in msg]
        if absent:
            check(False, "%s -> %s but message lacks %s: %s"
                  % (what, type(e).__name__, absent, msg[:220]))
        else:
            check(True, "%s -> %s; message: %s"
                  % (what, type(e).__name__, msg[:220]))
        return
    check(False, "%s -> nothing raised" % what)


def _with_ranks(table_name, value, fn):
    orig = dict(getattr(matrix, table_name))
    try:
        d = getattr(matrix, table_name)
        d.clear()
        d["fake.case"] = value
        return fn()
    finally:
        d = getattr(matrix, table_name)
        d.clear()
        d.update(orig)


def check_cell_cost_fortran_nonpositive_raises():
    raises(lambda: _with_ranks("FORTRAN_RANKS", 0,
                               lambda: run_e2e.cell_cost("fake.case", "fortran")),
           (ValueError,),
           "row112.i: cell_cost fortran with FORTRAN_RANKS['fake.case']=0",
           must_contain=("FORTRAN_RANKS", "fake.case"))
    raises(lambda: _with_ranks("FORTRAN_RANKS", -2,
                               lambda: run_e2e.cell_cost("fake.case", "fortran")),
           (ValueError,),
           "row112.i: cell_cost fortran with FORTRAN_RANKS['fake.case']=-2")
    got = _with_ranks("FORTRAN_RANKS", 4,
                      lambda: run_e2e.cell_cost("fake.case", "fortran"))
    check(got == 4, "row112.i: cell_cost fortran with a genuine positive "
                    "rank count (4) still returns it unchanged (got %r)" % got)


def check_cell_cost_jaxmpi_nonpositive_raises():
    raises(lambda: _with_ranks("PY_MPI_RANKS", 0,
                               lambda: run_e2e.cell_cost("fake.case",
                                                         "python-jax-mpi")),
           (ValueError,),
           "row112.i: cell_cost python-jax-mpi with PY_MPI_RANKS['fake.case']=0",
           must_contain=("PY_MPI_RANKS", "fake.case"))
    got = _with_ranks("PY_MPI_RANKS", 4,
                      lambda: run_e2e.cell_cost("fake.case", "python-jax-mpi"))
    check(got == 4, "row112.i: cell_cost python-jax-mpi with a genuine "
                    "positive rank count (4) still returns it (got %r)" % got)


def check_profile_ranks_fortran_nonpositive_raises():
    raises(lambda: _with_ranks("FORTRAN_RANKS", 0,
                               lambda: run_e2e.profile_ranks("fake.case",
                                                             "fortran")),
           (ValueError,),
           "row112.i: profile_ranks fortran with FORTRAN_RANKS['fake.case']=0",
           must_contain=("FORTRAN_RANKS", "fake.case"))
    got = _with_ranks("FORTRAN_RANKS", 8,
                      lambda: run_e2e.profile_ranks("fake.case", "fortran"))
    check(got == 8, "row112.i: profile_ranks fortran with a genuine "
                    "positive rank count (8) still returns it (got %r)" % got)


def check_profile_ranks_jaxmpi_nonpositive_raises():
    raises(lambda: _with_ranks("PY_MPI_RANKS", -1,
                               lambda: run_e2e.profile_ranks("fake.case",
                                                             "python-jax-mpi")),
           (ValueError,),
           "row112.i: profile_ranks python-jax-mpi with "
           "PY_MPI_RANKS['fake.case']=-1",
           must_contain=("PY_MPI_RANKS", "fake.case"))


def check_jax_measured_cores_ceiling_untouched():
    """The one site this fix deliberately leaves alone (mission item 112.i):
    JAX_MEASURED_CORES is a measured-cores figure ceil'd up, not a
    configured rank count -- cell_cost python-jax must still return a
    plain positive int via the untouched max(1, math.ceil(...))."""
    got = run_e2e.cell_cost("test.tpv8", "python-jax")
    check(isinstance(got, int) and got >= 1,
          "row112.i: cell_cost python-jax still returns a plain int >= 1 "
          "(got %r) -- this call path is untouched by the fix" % got)


CHECKS = [
    ("cell_cost-fortran", check_cell_cost_fortran_nonpositive_raises),
    ("cell_cost-jaxmpi", check_cell_cost_jaxmpi_nonpositive_raises),
    ("profile_ranks-fortran", check_profile_ranks_fortran_nonpositive_raises),
    ("profile_ranks-jaxmpi", check_profile_ranks_jaxmpi_nonpositive_raises),
    ("jax-ceiling-untouched", check_jax_measured_cores_ceiling_untouched),
]


def main():
    for tag, fn in CHECKS:
        print("\n-- %s --" % tag)
        fn()
    print()
    if FAILURES:
        print("FAIL test_perf_row112_e2e_ranks: %d check(s) failed" % len(FAILURES))
        return 1
    print("SUCCESS test_perf_row112_e2e_ranks: all %d checks passed" % len(CHECKS))
    return 0


if __name__ == "__main__":
    sys.exit(main())

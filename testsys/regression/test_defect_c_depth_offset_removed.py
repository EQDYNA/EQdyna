#! /usr/bin/env python3
"""
Defect C (board PR #164, 2026-10-09): setPlasticStress's depth argument used
to be `-0.5*(zline(iz)+zline(iz-1)) + 7.3215d0` in Fortran and the verbatim
mirror in Python -- an unexplained, unjustified additive constant that left a
nonzero sigma_zz at the true free surface in every Method-2 (C_elastic==0)
case (TPV13, TPV27, TPV30, drv.a6, drv.a6.v2), measured as a step-1 surface
residual acceleration of ~0.13-0.16 m/s^2 (EQDYNA_DUMP_EQUIL=1) that dropped
to ~0.00-0.05 m/s^2 after removal. Owner decision 2026-10-09, superseding
item 24(d) (2026-09-16, "do what fortran does ... not to be derived or
changed"): "it should be universal" / "we are setting stresses at the center
of each cell" -- depth is the element-centre depth everywhere, no offset, in
BOTH backends (rule: when you change physics, change it in both).

RULE 14a: this must come out BOTH ways.

  PASS case : the depth expression in both source files matches the
              fixed/expected one (no offset) -- the current, fixed tree.
  FAIL case : a synthetic copy of each source file with the historical
              "+ 7.3215" shift re-inserted is fed to the SAME checker and
              must be rejected.

Cheap (rule 9): grep/regex over already-committed source text, no solver run.
"""
import os
import re
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FORTRAN_SRC = os.path.join(REPO_ROOT, 'src', 'fortran', 'meshgen.f90')
PYTHON_SRC = os.path.join(REPO_ROOT, 'src', 'python', 'eqdyna', 'meshgen.py')

# The one call site in Fortran that feeds setPlasticStress its depth.
FORTRAN_CALL_RE = re.compile(
    r'call\s+setPlasticStress\(\s*([^,]+?)\s*,\s*elemCount\)')
# The two CODE sites that compute the same quantity (the vectorized
# build_elements assignment and the scalar mirror). The offset group is a
# GENERIC trailing "+ <number>" -- not literal "7.3215" -- so a regression
# that re-adds any other nonzero additive constant (not just the historical
# one) is caught too, not just a verbatim re-insertion of 7.3215.
PYTHON_DEPTH_RE = re.compile(
    r'depth(?:_val)?\s*=\s*(-0\.5\s*\*\s*\([^)]*\))(\s*\+\s*[0-9]+(?:\.[0-9]*)?)?')
# The two CODE sites (the vectorized build_elements assignment and the
# scalar mirror); the docstring's prose example (`-0.5*(zline[iz]+...`) is
# backtick-quoted text, not an assignment, and is NOT matched by
# PYTHON_DEPTH_RE at all (no `depth(_val) =` prefix there) -- it carries no
# live offset to catch a regression in, so it is intentionally left
# unchecked by this function rather than covered by some other check.
MIN_PYTHON_SITES = 2

EXPECTED_FORTRAN = '-0.5d0*(zline(iz)+zline(iz-1))'
FORBIDDEN_SUBSTR = '7.3215'


def check_fortran_text(text):
    """Return (ok, message). ok is False if the call site still carries the
    removed offset, or if the call site cannot be found at all (a change in
    shape this check did not anticipate must fail loudly, not pass silently)."""
    m = FORTRAN_CALL_RE.search(text)
    if m is None:
        return False, 'setPlasticStress call site not found at all'
    arg = m.group(1).strip()
    if FORBIDDEN_SUBSTR in arg:
        return False, 'setPlasticStress call site still carries the removed %s offset: %r' % (
            FORBIDDEN_SUBSTR, arg)
    if arg != EXPECTED_FORTRAN:
        return False, 'setPlasticStress call site does not match the expected depth expression: got %r, want %r' % (
            arg, EXPECTED_FORTRAN)
    return True, 'ok'


def check_python_text(text):
    sites = PYTHON_DEPTH_RE.findall(text)
    if len(sites) < MIN_PYTHON_SITES:
        return False, 'expected >= %d depth-expression CODE sites in meshgen.py, found %d' % (
            MIN_PYTHON_SITES, len(sites))
    for base, offset in sites:
        if offset and float(offset.replace('+', '').strip()) != 0.0:
            return False, (
                'a depth expression still carries a nonzero additive offset '
                '(not necessarily the historical 7.3215): %r%r' % (base, offset))
    return True, 'ok'


def main():
    failures = []

    with open(FORTRAN_SRC) as f:
        fortran_text = f.read()
    ok, msg = check_fortran_text(fortran_text)
    if not ok:
        failures.append('CURRENT TREE (Fortran, should PASS): %s' % msg)
    with open(PYTHON_SRC) as f:
        python_text = f.read()
    ok, msg = check_python_text(python_text)
    if not ok:
        failures.append('CURRENT TREE (Python, should PASS): %s' % msg)

    # Both-ways: re-insert the historical offset into each source TEXT (never
    # touching the files on disk) and confirm the SAME checker rejects it.
    reverted_fortran = fortran_text.replace(
        'call setPlasticStress(-0.5d0*(zline(iz)+zline(iz-1)), elemCount)',
        'call setPlasticStress(-0.5d0*(zline(iz)+zline(iz-1)) + 7.3215d0, elemCount)')
    if reverted_fortran == fortran_text:
        failures.append('could not construct the reverted-Fortran fixture -- '
                         'the call-site text this test substitutes against '
                         'has drifted out from under it')
    else:
        ok, msg = check_fortran_text(reverted_fortran)
        if ok:
            failures.append('REVERTED FIXTURE (Fortran, should FAIL): checker '
                             'passed a tree with the 7.3215 offset reinstated')

    reverted_python = python_text.replace(
        'depth = -0.5 * (zline[IZ] + zline[IZ - 1])',
        'depth = -0.5 * (zline[IZ] + zline[IZ - 1]) + 7.3215')
    if reverted_python == python_text:
        failures.append('could not construct the reverted-Python fixture -- '
                         'the depth-expression text this test substitutes '
                         'against has drifted out from under it')
    else:
        ok, msg = check_python_text(reverted_python)
        if ok:
            failures.append('REVERTED FIXTURE (Python, should FAIL): checker '
                             'passed a tree with the 7.3215 offset reinstated')

    # A SEPARATE both-ways fixture, a DIFFERENT additive constant (3.0, never
    # historically present) -- proves the checker catches ANY nonzero offset
    # reappearing, not merely a verbatim re-insertion of the literal 7.3215
    # this regression happened to be filed under.
    reverted_python_other = python_text.replace(
        'depth = -0.5 * (zline[IZ] + zline[IZ - 1])',
        'depth = -0.5 * (zline[IZ] + zline[IZ - 1]) + 3.0')
    if reverted_python_other == python_text:
        failures.append('could not construct the reverted-Python (other-constant) '
                         'fixture -- the depth-expression text this test substitutes '
                         'against has drifted out from under it')
    else:
        ok, msg = check_python_text(reverted_python_other)
        if ok:
            failures.append('REVERTED FIXTURE (Python, +3.0, should FAIL): checker '
                             'passed a tree with a NEW nonzero offset (not 7.3215) -- '
                             'the generalization this test is for did not land')

    if failures:
        for f in failures:
            print('FAIL: %s' % f)
        return 1
    print('PASS: Defect C offset absent in both backends (current tree), '
          'and the same checker rejects both backends\' reverted fixtures '
          '(rule 14a both-ways).')
    return 0


if __name__ == '__main__':
    sys.exit(main())

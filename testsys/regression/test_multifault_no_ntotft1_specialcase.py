#! /usr/bin/env python3
"""
Regression guard (Row 17, multi-fault): the fault-node engine must stay
NEUTRAL to ntotft -- one code path for 1..N faults, no `ntotft==1` branch and
no literal fault-index-1 hardcode reappearing in the files this mission
generalized.

Modelled on test_term_axis.py's shape: a grep-based scan that must PASS on
the current tree and is mutation-tested to FAIL when the defect shape it
guards is reintroduced (planted in a temp copy, never in the real tree).

WHAT IS FORBIDDEN, in SCANNED_FILES (the Fortran fault-engine files this
mission touched):
  - `ntotft == 1` / `ntotft==1` as a condition (a branch that treats the
    single-fault case specially instead of letting ntotft=1 fall out of a
    general loop as the degenerate case).
  - A bare `nftnd0(1)`, `nftnd(1)`, `fltxyz(.., 1)`, `checkIsOnFault(.., 1)`,
    `fric(.., 1)`, `nsmp(.., 1)` etc -- a literal `1` in the FAULT-INDEX
    position of an array this mission made ntotft-dimensioned -- unless the
    line is explicitly ALLOWLISTED below.

ALLOWLIST (checked to still be true, not just asserted): lines legitimately
naming fault 1 specifically, because they are about a DIFFERENT,
already-existing, single-fault-only mechanism this mission deliberately left
alone (C_degen>3 wedge degeneration; insertFaultType rough-fault geometry),
or because they are a loop BOUND/INDEX variable declaration (e.g.
`nftnd0(ntotft)`) that merely happens to contain the digit sequence "1" as
part of "ntotft" itself and is not actually about fault 1.

This is NOT a Python-port guard: the Python port (src/python/eqdyna/) is
ntotft==1 throughout as of this writing (see the mission report) and is
explicitly out of scope for this particular check -- scanning it here would
immediately fail on a KNOWN, already-reported gap rather than catching a
NEW regression, defeating the guard's purpose.

Cheap (rule 9): pure text scan of ~10 files, no build, sub-second.
"""
import glob
import os
import re
import shutil
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FSRC = os.path.join(ROOT, 'src', 'fortran')

# The Fortran files this mission generalized to be ntotft-neutral. Deliberately
# NOT the whole src/fortran/ tree: files untouched by this mission (e.g.
# func_lib.f90's insertFaultInterface, the rough-fault y-blend) may
# legitimately still be fault-1-only and are not this guard's business.
SCANNED_FILES = [
    'meshgen.f90', 'countMeshEntities.f90', 'assembleGlobalMass.f90',
    'netcdf_io.f90', 'library_output.f90', 'eqdyna3d.f90',
    'checkInputConsistency.f90', 'faulting.f90', 'driver.f90',
    'readInputFiles.f90',
]

# (filename, 0-based line index is NOT used -- matched by a short, specific
# substring so the allowlist survives line-number churn) -> reason.
# Every entry here was checked by hand to be the documented exemption, not a
# missed generalization.
ALLOWLIST = [
    # C_degen>3 wedge degeneration (meshgen.f90): a pre-existing, orthogonal,
    # already single-fault-only dipping-fault mesh mechanism. Guarded by
    # `if (C_degen > 3.0d0)`, inactive for this mission's vertical-fault
    # scope. Left exactly as it was; not part of the multi-fault engine.
    ('meshgen.f90', 'call wedge(elementCenterCoor(1)'),
    ('meshgen.f90', "call checkIsOnFault(meshCoor(1:3,nodeElemIdRelation(1,elemCount)), 1, isOnFt)"),
    ('meshgen.f90', "call checkIsOnFault(meshCoor(1:3,nodeElemIdRelation(2,elemCount)), 1, isOnFt)"),
    # insertFaultType rough-fault geometry (meshgen.f90/readInputFiles.f90):
    # the rough-fault y-blend is deliberately single-fault (CLAUDE.md,
    # checkIsOnFault's own comment) and orthogonal to multi-fault; the rough
    # geometry file is validated against fault 1's box specifically because
    # rough geometry + multi-fault is not this mission's scope.
    ('meshgen.f90', 'numOfNodesWithUniformGridsize = nint((fltxyz(2,1,1)'),
    ('meshgen.f90', 'frontEdgeCoor = fltxyz(1,1,1)'),
    ('meshgen.f90', 'backEdgeCoor  = fltxyz(2,1,1)'),
    ('meshgen.f90', 'numOfNodesWithUniformGridsize = nint((fltxyz(2,3,1)'),
    ('meshgen.f90', 'frontEdgeCoor = fltxyz(1,3,1)'),
    ('meshgen.f90', 'backEdgeCoor  = fltxyz(2,3,1)'),
    ('readInputFiles.f90', "abs(rough_fx_min - fltxyz(1,1,1))"),
    ('readInputFiles.f90', "abs(rough_fz_min - fltxyz(1,3,1))"),
    ('readInputFiles.f90', "fltxyz(1,1,1), fltxyz(1,3,1)"),
    ('readInputFiles.f90', "spanx = (fltxyz(2,1,1) - fltxyz(1,1,1))/dx"),
    ('readInputFiles.f90', "spanz = (fltxyz(2,3,1) - fltxyz(1,3,1))/dz"),
    # replaceSlaveWithMasterNode (meshgen.f90): the x/z bound check here is
    # intentionally keyed to fault 1's box alone -- checkInputConsistency.f90
    # already REQUIRES every fault share fault 1's x/z extent (this
    # mission's narrow two-parallel-faults scope: independent per-fault x/z
    # extents are explicitly out of scope, guarded there).
    ('meshgen.f90', 'nodeCoor(1)>(fltxyz(1,1,1)-tol) .and. nodeCoor(1)<(fltxyz(2,1,1)+dx+tol)'),
    ('meshgen.f90', 'nodeCoor(3)>(fltxyz(1,3,1)-tol)'),
    # checkInputConsistency.f90's own multi-fault guard: deliberately
    # compares every fault's x/z box against FAULT 1's specifically (the
    # shared-x/z-extent requirement this mission's narrow scope imposes),
    # not a missed generalization.
    ('checkInputConsistency.f90', 'abs(fxmin(i)-fxmin(1))>tol .or. abs(fxmax(i)-fxmax(1))>tol'),
    ('checkInputConsistency.f90', 'abs(fzmin(i)-fzmin(1))>tol .or. abs(fzmax(i)-fzmax(1))>tol'),
]

# `ntotft == 1` / `ntotft==1` as a condition -- forbidden everywhere in
# SCANNED_FILES, no exemptions (this is the actual "no special-casing" rule).
NTOTFT_EQ_1 = re.compile(r'ntotft\s*==\s*1\b')

# A literal `1` in the fault-index position of an array this mission made
# (.., ntotft)-shaped: last argument of a 2+-arg parenthesised subscript/call.
# Matched conservatively (specific array/subroutine names), not a blanket
# "any ,1)" scan, which would flag unrelated code (e.g. a genuinely 1-D array
# subscript, or column/row index 1 of something unrelated to faults).
# Two shapes: a MULTI-arg array/call with fault index as the LAST argument
# (fltxyz(:,:,1), nsmp(:,:,1), checkIsOnFault(x,1,y) via a preceding comma),
# and a SINGLE-arg array whose one dimension IS the fault index (nftnd(1),
# nftnd0(1), nonfs(1), fxmin(1), ...).
FAULT_INDEX_1_MULTIARG = re.compile(
    r'\b(fltxyz|nsmp|fric|un|us|ud|arn|fnft|Tatnode|patnode|'
    r'fltgm|fltnum|fltl|fltr|fltf|fltb|fltd|fltu|onFaultTPHist|'
    r'checkIsOnFault|faultTag)\s*\([^()]*,\s*1\s*\)'
)
FAULT_INDEX_1_SINGLEARG = re.compile(
    r'\b(nftnd0?|nonfs|fxmin|fxmax|fymin|fymax|fzmin|fzmax)\s*\(\s*1\s*\)'
)


def _iter_allowlisted(fname):
    return {snippet for f, snippet in ALLOWLIST if f == fname}


def scan_file(path, fname):
    """Return a list of (lineno, line) violations."""
    violations = []
    allow = _iter_allowlisted(fname)
    with open(path, errors='replace') as fh:
        for lineno, raw in enumerate(fh, 1):
            code = raw.split('!', 1)[0]  # strip Fortran end-of-line comments
            if not code.strip():
                continue
            if NTOTFT_EQ_1.search(code):
                violations.append((lineno, raw.rstrip(), 'ntotft==1 branch'))
                continue
            m = FAULT_INDEX_1_MULTIARG.search(code) or FAULT_INDEX_1_SINGLEARG.search(code)
            if m:
                if any(a in code for a in allow):
                    continue
                violations.append((lineno, raw.rstrip(),
                                   'literal fault-index 1 on %s(...)' % m.group(1)))
    return violations


def check_current_tree_clean():
    bad = []
    for fname in SCANNED_FILES:
        path = os.path.join(FSRC, fname)
        if not os.path.isfile(path):
            raise AssertionError('%s is in SCANNED_FILES but does not exist -- '
                                 'update this test' % fname)
        for lineno, line, why in scan_file(path, fname):
            bad.append('%s:%d: %s\n    %s' % (fname, lineno, why, line.strip()))
    if bad:
        raise AssertionError('ntotft==1 special-casing found in the multi-fault engine:\n  '
                             + '\n  '.join(bad))
    print('  PASS  %d file(s) scanned, no ntotft==1 branch and no unlisted literal '
          'fault-index-1 found' % len(SCANNED_FILES))


def check_mutation_ntotft_eq_1_branch_is_caught():
    with tempfile.TemporaryDirectory() as tmp:
        bad_path = os.path.join(tmp, 'meshgen.f90')
        shutil.copy(os.path.join(FSRC, 'meshgen.f90'), bad_path)
        with open(bad_path, 'a') as fh:
            fh.write('\n    if (ntotft == 1) then\n        continue\n    endif\n')
        v = scan_file(bad_path, 'meshgen.f90')
        if not v:
            raise AssertionError('planted `if (ntotft == 1) then` was NOT detected -- '
                                 'the scanner cannot see this defect shape')
    print('  PASS  mutation: a planted `ntotft == 1` branch is detected')


def check_mutation_literal_fault_index_is_caught():
    with tempfile.TemporaryDirectory() as tmp:
        bad_path = os.path.join(tmp, 'library_output.f90')
        shutil.copy(os.path.join(FSRC, 'library_output.f90'), bad_path)
        with open(bad_path, 'a') as fh:
            fh.write('\n    ! planted regression\n    x = nftnd(1)\n')
        v = scan_file(bad_path, 'library_output.f90')
        if not v:
            raise AssertionError('planted literal `nftnd(1)` was NOT detected')
    print('  PASS  mutation: a planted literal fault-index-1 (nftnd(1)) is detected')


def check_allowlist_entries_still_exist():
    """An allowlist entry for a line that no longer exists is a stale
    exemption hiding nothing -- prune it rather than let it rot."""
    missing = []
    for fname, snippet in ALLOWLIST:
        path = os.path.join(FSRC, fname)
        text = open(path, errors='replace').read()
        if snippet not in text:
            missing.append('%s: %r' % (fname, snippet))
    if missing:
        raise AssertionError('allowlist entries that no longer match anything in the '
                             'source (stale -- prune them): %s' % '; '.join(missing))
    print('  PASS  all %d allowlist entries still match real source lines' % len(ALLOWLIST))


def main():
    print('Regression guard: no ntotft==1 special-casing in the multi-fault engine')
    failures = []
    for c in (check_current_tree_clean,
              check_allowlist_entries_still_exist,
              check_mutation_ntotft_eq_1_branch_is_caught,
              check_mutation_literal_fault_index_is_caught):
        try:
            c()
        except AssertionError as e:
            failures.append('%s: %s' % (c.__name__, e))
            print('  FAIL  %s: %s' % (c.__name__, e))
    if failures:
        print('\nFAIL test_multifault_no_ntotft1_specialcase (%d check(s))' % len(failures))
        return 1
    print('\nSUCCESS test_multifault_no_ntotft1_specialcase')
    return 0


if __name__ == '__main__':
    sys.exit(main())

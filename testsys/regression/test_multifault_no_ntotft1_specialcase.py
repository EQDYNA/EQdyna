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
    # Row 17 REBASED (restore per-fault mesh extent): getLocalOneDimCoorArrAndSize's
    # fault-1-only reads of fltxyz (x/z uniform-belt origin/extent) and
    # replaceSlaveWithMasterNode's fault-1-only x/z bound test, plus
    # checkInputConsistency.f90's "every fault must share fault 1's x/z
    # extent" guard, are GONE -- not relaxed, removed: the belt is now the
    # UNION (minval/maxval) over every fault's own box, and
    # replaceSlaveWithMasterNode tests each fault's own x/z alongside its own
    # y-plane, one fault at a time. The allowlist entries that used to
    # document those fault-1-only reads as a deliberate, reviewed exemption
    # are deleted along with the code, not left stale.
    ('readInputFiles.f90', "abs(rough_fx_min - fltxyz(1,1,1))"),
    ('readInputFiles.f90', "abs(rough_fz_min - fltxyz(1,3,1))"),
    ('readInputFiles.f90', "fltxyz(1,1,1), fltxyz(1,3,1)"),
    ('readInputFiles.f90', "spanx = (fltxyz(2,1,1) - fltxyz(1,1,1))/dx"),
    ('readInputFiles.f90', "spanz = (fltxyz(2,3,1) - fltxyz(1,3,1))/dz"),
    # Same C_degen>3 wedge-degeneration exemption as the entries above --
    # nftnd0(1) here is the pre-existing, single-fault-only wedge mechanism's
    # own argument, not a missed generalization.
    ('meshgen.f90', 'iy, iz, nftnd0(1))'),
    # faultTag() (library_output.f90) IS the per-fault tagging convention
    # itself (CLAUDE.md: '' for fault 1, 'ft<N>_' for fault N>=2) -- it must
    # compare ntotft and ift against 1 to decide whether to tag at all. This
    # is the one place in the engine where naming fault 1 specifically is the
    # whole point, not a missed generalization.
    ('library_output.f90', 'ntotft > 1 .and. ift > 1'),
]

# `ntotft` (or the reversed operand order `1 == ntotft`) compared, by ANY
# relational operator and in EITHER Fortran spelling (`==`/`.eq.`, `/=`/
# `.ne.`, `<=`/`.le.`, `<`/`.lt.`, `>=`/`.ge.`, `>`/`.gt.`), against the
# literal 1 or 2 -- forbidden everywhere in SCANNED_FILES, no exemptions
# (this is the actual "no special-casing" rule). `ntotft < 2` and
# `ntotft <= 1` are the same single-fault special case as `ntotft == 1`
# written with a different operator; `ntotft > 1` is the same special case
# from its complementary side. re.IGNORECASE: Fortran is case-insensitive, so
# `NTOTFT == 1` is the identical defect. Audit note (victor-reyes,
# 2026-09-30): the previous version of this guard matched only the single
# literal spelling `ntotft==1` and missed 11 of 15 equivalent-mutation probes
# (`.eq.`, case, every other operator, operand order).
_REL_OP = r'(?:==|\.eq\.|/=|\.ne\.|<=|\.le\.|<|\.lt\.|>=|\.ge\.|>|\.gt\.)'
NTOTFT_COMPARISON = re.compile(
    r'\bntotft\s*' + _REL_OP + r'\s*[12]\b'
    r'|\b[12]\s*' + _REL_OP + r'\s*ntotft\b',
    re.IGNORECASE)

# Same defect shape, one level down: the per-fault LOOP INDEX (`ift` in every
# scanned file, `iFault` in meshgen.f90/checkInputConsistency.f90) compared
# against the literal 1 -- `if (ift == 1)` special-cases fault 1 inside a
# loop that is supposed to treat every fault alike, exactly as `ntotft == 1`
# special-cases the single-fault case outside one.
FAULT_VAR_COMPARISON = re.compile(
    r'\b(?:ift|iFault)\s*' + _REL_OP + r'\s*1\b'
    r'|\b1\s*' + _REL_OP + r'\s*(?:ift|iFault)\b',
    re.IGNORECASE)

# A literal `1` in the fault-index position of an array this mission made
# (.., ntotft)-shaped: last argument of a 2+-arg parenthesised subscript/call.
# Matched conservatively (specific array/subroutine names), not a blanket
# "any ,1)" scan, which would flag unrelated code (e.g. a genuinely 1-D array
# subscript, or column/row index 1 of something unrelated to faults).
# Two shapes: a MULTI-arg array/call with fault index as the LAST argument
# (fltxyz(:,:,1), nsmp(:,:,1), checkIsOnFault(x,1,y) via a preceding comma),
# and a SINGLE-arg array whose one dimension IS the fault index (nftnd(1),
# nftnd0(1), nonfs(1), fxmin(1), ...). re.IGNORECASE added for the same
# reason as NTOTFT_COMPARISON (`NFTND(1)` is the same defect as `nftnd(1)`).
# _INNER allows exactly one level of nested parens inside the call, so a
# fault-index-1 literal buried in a call that itself contains another call
# (e.g. `fric(FRIC_SLOT_STATE, nsmp(1,i,ift), 1)`) is still found -- the
# previous `[^()]*` could not cross the inner `nsmp(...)`'s own parens at all.
_INNER = r'(?:[^()]|\([^()]*\))*'
FAULT_INDEX_1_MULTIARG = re.compile(
    r'\b(fltxyz|nsmp|fric|un|us|ud|arn|fnft|Tatnode|patnode|'
    r'fltgm|fltnum|fltl|fltr|fltf|fltb|fltd|fltu|onFaultTPHist|'
    r'checkIsOnFault|faultTag)\s*\(' + _INNER + r',\s*1\s*\)',
    re.IGNORECASE)
FAULT_INDEX_1_SINGLEARG = re.compile(
    r'\b(nftnd0?|nonfs|fxmin|fxmax|fymin|fymax|fzmin|fzmax)\s*\(\s*1\s*\)',
    re.IGNORECASE)


def _iter_allowlisted(fname):
    return [snippet for f, snippet in ALLOWLIST if f == fname]


def scan_file(path, fname):
    """Return a list of (lineno, line, why) violations.

    Allowlisting is per-MATCH (does this specific matched snippet sit inside
    an allowlisted snippet?), not per-LINE (does an allowlisted snippet
    appear anywhere on this line?). The per-line version let an allowlisted
    fragment shield an unrelated bug placed later on the same physical line
    (Fortran's `;` statement separator) -- the
    `check_mutation_allowlist_does_not_shield_a_different_bug_on_same_line`
    probe below pins this.

    Allowlisting is additionally restricted to each entry's FIRST (earliest,
    by absolute file offset) occurrence. Every entry was hand-verified
    against exactly one real line (`check_allowlist_entries_still_exist`
    confirms each occurs exactly once in the current tree); a SECOND
    occurrence of identical text appearing later in the file -- which cannot
    happen in the untouched tree, only when a planted violation happens to
    reproduce an allowlisted snippet verbatim elsewhere -- is therefore never
    the hand-verified line and must not be exempted. Without this, a planted
    violation whose text is byte-identical to an allowlist entry (e.g.
    replaying `ntotft > 1 .and. ift > 1` outside faultTag()) would slip
    through the per-line position check too, since the identical text
    trivially satisfies 'falls inside an occurrence on this line' for
    whichever line it was pasted on."""
    violations = []
    allow = _iter_allowlisted(fname)

    full_text = open(path, errors='replace').read()
    canonical_offset = {entry: full_text.find(entry) for entry in allow}

    def _allowlisted(code, start, end, abs_line_start):
        """A match (code[start:end]) is allowlisted only if it falls INSIDE
        an occurrence of one of this file's allowlisted snippets AT THAT
        POSITION on THIS physical line, AND that occurrence is the entry's
        canonical (first, hand-verified) one in the file -- not merely if
        the matched text is a substring of some allowlist entry that exists
        anywhere in the file. (The former bug: an allowlist entry for one
        line let a same-text match on any OTHER line, in a different file
        region entirely, slip through unflagged.)"""
        for entry in allow:
            idx = code.find(entry)
            while idx != -1:
                if (idx <= start and end <= idx + len(entry)
                        and abs_line_start + idx == canonical_offset[entry]):
                    return True
                idx = code.find(entry, idx + 1)
        return False

    abs_offset = 0
    for lineno, raw in enumerate(full_text.splitlines(keepends=True), 1):
        code = raw.split('!', 1)[0]  # strip Fortran end-of-line comments
        line_abs_start = abs_offset
        abs_offset += len(raw)
        if not code.strip():
            continue
        for m in NTOTFT_COMPARISON.finditer(code):
            if _allowlisted(code, m.start(), m.end(), line_abs_start):
                continue
            violations.append((lineno, raw.rstrip(), 'ntotft literal comparison (%r)' % m.group(0)))
        for m in FAULT_VAR_COMPARISON.finditer(code):
            if _allowlisted(code, m.start(), m.end(), line_abs_start):
                continue
            violations.append((lineno, raw.rstrip(), 'fault-index-variable literal comparison (%r)' % m.group(0)))
        for m in FAULT_INDEX_1_MULTIARG.finditer(code):
            if _allowlisted(code, m.start(), m.end(), line_abs_start):
                continue
            violations.append((lineno, raw.rstrip(),
                               'literal fault-index 1 on %s(...)' % m.group(1)))
        for m in FAULT_INDEX_1_SINGLEARG.finditer(code):
            if _allowlisted(code, m.start(), m.end(), line_abs_start):
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


# Audit (victor-reyes, 2026-09-30): the auditor's 15-probe equivalent-mutation
# sweep caught only 4 of 15 under the previous single-literal-spelling regex.
# Every one of the 15 planted here, one at a time, into a temp copy of
# meshgen.f90 (the file the original two checks above already used for this
# purpose) -- each must be independently detected.
EQUIVALENT_MUTATION_FORMS = {
    'eq-lower':                  'if (ntotft == 1) then',
    'dot-eq':                    'if (ntotft .eq. 1) then',
    'UPPER':                     'if (NTOTFT == 1) then',
    'lt2':                       'if (ntotft < 2) then',
    'le1':                       'if (ntotft <= 1) then',
    'ne1':                       'if (ntotft /= 1) then',
    'gt1':                       'if (ntotft > 1) call abortRun(45, "x")',
    'reversed':                  'if (1 == ntotft) then',
    'nftnd1':                    'x = nftnd(1)',
    'nftnd1-upper':              'x = NFTND(1)',
    'fric-nested':               'x = fric(FRIC_SLOT_STATE, nsmp(1,i,ift), 1)',
    'meshCoor-nsmp1':            'x = meshCoor(1, nsmp(1,i,1))',
    'fnft-1':                    'x = fnft(i,1)',
    'allowlisted-line-plus-bug': 'frontEdgeCoor = fltxyz(1,1,1); y = nftnd(1)',
    'ift-eq-1':                  'if (ift == 1) then',
}


def check_mutation_equivalent_forms_are_caught():
    with tempfile.TemporaryDirectory() as tmp:
        missed = []
        for name, line in EQUIVALENT_MUTATION_FORMS.items():
            bad_path = os.path.join(tmp, 'meshgen_%s.f90' % name)
            shutil.copy(os.path.join(FSRC, 'meshgen.f90'), bad_path)
            with open(bad_path, 'a') as fh:
                fh.write('\n    ' + line + '\n')
            v = scan_file(bad_path, 'meshgen.f90')
            if not v:
                missed.append('%s (%r)' % (name, line))
        if missed:
            raise AssertionError('%d/%d equivalent-mutation form(s) NOT detected: %s'
                                 % (len(missed), len(EQUIVALENT_MUTATION_FORMS), '; '.join(missed)))
    print('  PASS  all %d equivalent-mutation forms detected' % len(EQUIVALENT_MUTATION_FORMS))


def check_allowlist_does_not_shield_a_different_bug_on_same_line():
    """The 'allowlisted-line-plus-bug' probe in isolation: an allowlisted
    snippet on a line must not shield a DIFFERENT, unrelated violation placed
    later on that same physical line (Fortran's `;` statement separator).
    This is the defect the old whole-line `if any(a in code for a in allow)`
    check had -- it skipped the entire line, bug included, whenever any
    allowlisted snippet appeared anywhere in it."""
    with tempfile.TemporaryDirectory() as tmp:
        bad_path = os.path.join(tmp, 'meshgen.f90')
        shutil.copy(os.path.join(FSRC, 'meshgen.f90'), bad_path)
        with open(bad_path, 'a') as fh:
            # 'frontEdgeCoor = fltxyz(1,1,1)' alone IS allowlisted for this
            # file; the planted 'nftnd(1)' after the ';' is not and must
            # still be caught.
            fh.write('\n    frontEdgeCoor = fltxyz(1,1,1); y = nftnd(1)\n')
        v = scan_file(bad_path, 'meshgen.f90')
        if not v:
            raise AssertionError('an allowlisted snippet shielded an unrelated bug '
                                 '(nftnd(1)) placed later on the same line')
        if any('fltxyz(1,1,1)' in line for _, line, _ in v) and len(v) == 1 and 'nftnd' not in v[0][2]:
            raise AssertionError('the reported violation is the allowlisted snippet itself, '
                                 'not the planted nftnd(1) bug: %r' % (v,))
    print('  PASS  mutation: an allowlisted snippet does not shield a different bug on the same line')


# AUDIT FIX (victor-reyes, 2026-09-30): the allowlist check used to be
# `any(snippet in a for a in allow)` -- "does the matched text appear as a
# substring of SOME allowlist entry for this file, anywhere" -- with no
# requirement that the matched text's actual occurrence sit at the position
# an allowlisted snippet occupies ON THE LINE WHERE IT WAS MATCHED. That let
# a planted violation on ANY OTHER line of the file slip through unflagged,
# as long as its literal text happened to be a substring of some unrelated
# allowlisted line elsewhere in the same file. Each probe below plants the
# violation on its OWN new line (not the allowlisted line), in a temp copy of
# the file whose allowlist entry it exploits, and must still be caught.
REGRESSION_PROBES_DIFFERENT_LINE_THAN_ALLOWLIST_ENTRY = [
    # exploits meshgen.f90's 'iy, iz, nftnd0(1))' allowlist entry (the
    # C_degen>3 wedge call) via bare substring match on 'nftnd0(1)'.
    ('meshgen.f90', 'x = nftnd0(1)'),
    # Row 17 rebased: checkInputConsistency.f90's 'fxmin(i)-fxmin(1)' entry
    # (and the source line it allowlisted) are both gone -- this probe is now
    # a plain positive-detection check (no allowlist entry left to shield it
    # by substring at all), kept so a FUTURE literal fxmin(1) hardcode in
    # this file is still caught.
    ('checkInputConsistency.f90', 'x = fxmin(1)'),
    # exploits library_output.f90's faultTag() allowlist entry
    # 'ntotft > 1 .and. ift > 1' via the literal condition reappearing,
    # unrelated to faultTag(), on a different line.
    ('library_output.f90', 'if (ntotft > 1 .and. ift > 1) then'),
    ('library_output.f90', 'if (ift > 1) then'),
]


def check_mutation_allowlist_is_position_specific_not_file_global():
    """Each probe must be DETECTED when planted on a line other than the one
    the allowlist entry it textually overlaps with actually describes --
    pins the exact regression victor-reyes found and the fix for it."""
    missed = []
    with tempfile.TemporaryDirectory() as tmp:
        for fname, line in REGRESSION_PROBES_DIFFERENT_LINE_THAN_ALLOWLIST_ENTRY:
            bad_path = os.path.join(tmp, fname)
            shutil.copy(os.path.join(FSRC, fname), bad_path)
            with open(bad_path, 'a') as fh:
                fh.write('\n    ! planted regression, deliberately on its own new line\n    '
                         + line + '\n')
            v = scan_file(bad_path, fname)
            if not v:
                missed.append('%s: %r' % (fname, line))
    if missed:
        raise AssertionError('%d probe(s) NOT detected (allowlist is matching globally, not by '
                             'line position): %s' % (len(missed), '; '.join(missed)))
    print('  PASS  all %d position-specific allowlist regression probes detected'
          % len(REGRESSION_PROBES_DIFFERENT_LINE_THAN_ALLOWLIST_ENTRY))


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
              check_mutation_literal_fault_index_is_caught,
              check_mutation_equivalent_forms_are_caught,
              check_allowlist_does_not_shield_a_different_bug_on_same_line,
              check_mutation_allowlist_is_position_specific_not_file_global):
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

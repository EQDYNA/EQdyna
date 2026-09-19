#! /usr/bin/env python3
"""
Regression guard: CI's matrix jobs must cover matrix.CI_CELLS exactly --
no cell dropped, none silently run twice (rules 2, 6).

The e2e-ci selection is split across several jobs' `--backends`/`--cases`
filters so each gets its own runner. Splitting a declared cell list across
job command lines is exactly the shape of defect that drifts silently: a
later edit to matrix.CI_CELLS has no mechanical reason to also update the
workflow's filters, and vice versa -- a dropped cell would just quietly stop
being gated with a green run to show for it.

WHAT THIS PINS. It parses test.yml for every `run_e2e.py --ci ...` invocation,
applies the SAME --backends/--cases filtering run_e2e.py's own select() uses,
and asserts:
  1. at least one such invocation exists (a workflow that stopped using --ci
     entirely would otherwise pass this file vacuously);
  2. the union of what all of them select is EXACTLY matrix.CI_CELLS -- no
     more (that would be a cell CI silently started grading that the table
     never declared), no less (a cell the table declares but no job runs);
  3. no cell is selected by more than one invocation (redundant work is the
     opposite of what this split exists for, and a duplicate loosely masks a
     missing exclusive assignment).

Cheap (rule 9): text parsing plus a set comparison against an already-imported
table, well under 1 s. Exits non-zero on any failure.
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
WORKFLOW = os.path.join(ROOT, '.github', 'workflows', 'test.yml')

if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from testsys import matrix  # noqa: E402


def find_ci_invocations(text):
    """Every line invoking `run_e2e.py ... --ci ...`, with its --backends and
    --cases values (None if the flag is absent -- meaning 'all').

    Strips everything from the first '#' onward before matching: this file's
    own explanatory comments describe `run_e2e.py --ci ...` invocations in
    prose right next to the real ones, and matching prose instead of code is
    the exact false-positive failure mode
    testsys/regression/test_sweep_core_budget.py's own comment-stripping
    already warns about for a different guard. A line that is ALL comment
    (or blank after stripping) simply has nothing left to match.
    """
    invocations = []
    for raw in text.splitlines():
        line = raw.split('#', 1)[0]
        if 'run_e2e.py' not in line or '--ci' not in line:
            continue
        m_backends = re.search(r'--backends\s+([\w,\-]+)', line)
        m_cases = re.search(r'--cases\s+([\w,.]+)', line)
        backends = m_backends.group(1).split(',') if m_backends else None
        cases = m_cases.group(1).split(',') if m_cases else None
        invocations.append((line.strip(), cases, backends))
    return invocations


def filtered_cells(cases, backends):
    wanted = list(matrix.CI_CELLS)
    if cases:
        want_cases = set(cases)
        wanted = [c for c in wanted if c[0] in want_cases]
    if backends:
        want_backends = set(backends)
        wanted = [c for c in wanted if c[1] in want_backends]
    return set(wanted)


def main():
    print('Regression guard: CI matrix jobs must cover matrix.CI_CELLS exactly')
    text = open(WORKFLOW, errors='replace').read()
    invocations = find_ci_invocations(text)
    if not invocations:
        print('FAIL: no `run_e2e.py --ci` invocation found in %s -- either '
              'CI stopped gating the e2e-ci selection entirely, or it was '
              'rewritten to not use --ci, and this guard cannot see what it '
              'covers' % WORKFLOW)
        return 1

    print('  found %d invocation(s):' % len(invocations))
    all_cells = set(matrix.CI_CELLS)
    seen = {}
    problems = []
    for line, cases, backends in invocations:
        cells = filtered_cells(cases, backends)
        print('    %s -> %d cell(s)' % (line, len(cells)))
        for cell in cells:
            seen.setdefault(cell, []).append(line)

    covered = set(seen)
    missing = all_cells - covered
    extra = covered - all_cells
    duplicated = {cell: lines for cell, lines in seen.items() if len(lines) > 1}

    if missing:
        problems.append(
            'matrix.CI_CELLS has %d cell(s) no workflow invocation selects: %s'
            % (len(missing), ', '.join('%s x %s' % c for c in sorted(missing))))
    if extra:
        problems.append(
            'workflow invocation(s) select %d cell(s) NOT in matrix.CI_CELLS: '
            '%s' % (len(extra), ', '.join('%s x %s' % c for c in sorted(extra))))
    if duplicated:
        for cell, lines in duplicated.items():
            problems.append(
                '%s x %s is selected by %d invocations (wasted runner time, '
                'and masks which job is authoritative for it): %r'
                % (cell[0], cell[1], len(lines), lines))

    if problems:
        print('\nFAIL: %d problem(s):' % len(problems))
        for p in problems:
            print('  - %s' % p)
        return 1
    print('\nPASS: workflow invocations cover matrix.CI_CELLS exactly, no '
          'overlap (%d cells total)' % len(all_cells))
    return 0


if __name__ == '__main__':
    sys.exit(main())

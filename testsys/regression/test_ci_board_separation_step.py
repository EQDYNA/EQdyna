#! /usr/bin/env python3
"""
Regression guard: CI must still RUN the board/code separation check, in a job
that can actually resolve a commit range (rule 21c, pathway items 75/80).

THE INCIDENT SHAPE. `testsys/hooks/pre-commit` cannot see a `git commit
--amend`, so the merge-boundary half of rule 21c lives in a workflow step --
the `Check board/code commit separation (rule 21c)` step of test.yml's
`unit-regression` job, which runs `python3 testsys/check_board_separation.py
<range>`. Every OTHER mechanism that landed with it is guarded by a
regression test; this one was guarded by nothing. Three edits would disarm it
and leave every tier exiting 0:

  1. delete the step (or the `check_board_separation.py` invocation inside it);
  2. put `fetch-depth: 1` (or drop `fetch-depth` entirely, whose DEFAULT is 1)
     on that job's checkout -- the range endpoints then do not exist in the
     clone, and the step either dies on an unreachable sha or falls back to a
     range it cannot resolve;
  3. delete `testsys/check_board_separation.py`, or drop its exec bit, so the
     step invokes nothing.

WHAT THIS PINS, in order:
  1. exactly one step in the whole workflow invokes check_board_separation.py,
     with a range ARGUMENT, and that invocation is real YAML content -- this
     file PARSES test.yml structurally (jobs -> steps -> `run:` block scalars)
     rather than grepping the raw text, because test.yml's own comments
     describe this step in prose two ways ("the board-separation step below",
     "python3 testsys/check_board_separation.py") and a grep would be
     satisfied by the comments alone after the step was deleted;
  2. the step carries no `if:` -- `if: false` disarms a step while leaving it
     visible and the job green;
  3. the job containing it checks out with `fetch-depth: 0`;
  4. testsys/check_board_separation.py exists, is executable, and still
     REFUSES a missing range with exit 2 (its documented contract) -- proof
     the step's target is runnable, not just present.

Why a sibling of test_ci_workflow_coverage.py rather than an extension of it:
that file answers one question (does the union of `run_e2e.py --ci`
invocations equal matrix.CI_CELLS) with one mechanism (line regex), and this
one needs structural parsing plus a subprocess probe of a different artifact.
One guard per incident (rule 10) keeps the FAIL line unambiguous.

Cheap (rule 9): one file parse plus one sub-second subprocess. Exits non-zero
on any failure (rule 2).
"""
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
WORKFLOW = os.path.join(ROOT, '.github', 'workflows', 'test.yml')
CHECKER = os.path.join(ROOT, 'testsys', 'check_board_separation.py')
CHECKER_REL = os.path.join('testsys', 'check_board_separation.py')

KEY_RE = re.compile(r'^([A-Za-z0-9_.\-]+):(?:\s+(.*))?$')


class YamlSubsetError(Exception):
    """test.yml grew a construct this parser does not model. Raised, never
    swallowed: a parser that silently returns a partial document would let
    this guard pass on a file it did not actually understand."""


# ---------------------------------------------------------------- parser ---
# A deliberately small, strict indentation parser for the subset test.yml
# uses: block mappings, block sequences, block scalars (`|`), `#` comments.
# Stdlib only -- PyYAML is not in CI's `pip install` line, and adding it there
# would be a change to the workflow this file exists to guard.

def _indent(line):
    return len(line) - len(line.lstrip(' '))


def _significant(lines, i):
    while i < len(lines) and (not lines[i].strip() or lines[i].lstrip().startswith('#')):
        i += 1
    return i


def _strip_trailing_comment(value):
    # Only ` #` (space-hash) starts a YAML trailing comment; a bare '#' inside
    # a token does not. No quoted '#' occurs in this file's scalars.
    cut = value.find(' #')
    value = (value[:cut] if cut >= 0 else value).strip()
    # Unquote, so `fetch-depth: "0"` and `fetch-depth: 0` compare the same.
    if len(value) >= 2 and value[0] == value[-1] and value[0] in ('"', "'"):
        value = value[1:-1]
    return value


def _read_block_scalar(lines, i, parent_indent):
    body = []
    while i < len(lines):
        if not lines[i].strip():
            body.append('')
            i += 1
            continue
        if _indent(lines[i]) <= parent_indent:
            break
        body.append(lines[i])
        i += 1
    while body and not body[-1].strip():
        body.pop()
    widths = [_indent(b) for b in body if b.strip()]
    pad = min(widths) if widths else 0
    return '\n'.join(b[pad:] if b.strip() else '' for b in body), i


def _parse_block(lines, i, min_indent):
    i = _significant(lines, i)
    if i >= len(lines) or _indent(lines[i]) < min_indent:
        return None, i
    cur = _indent(lines[i])
    if lines[i].lstrip().startswith('- '):
        seq = []
        while True:
            i = _significant(lines, i)
            if (i >= len(lines) or _indent(lines[i]) != cur
                    or not lines[i].lstrip().startswith('- ')):
                break
            rest = lines[i].lstrip()[2:]
            item_indent = cur + 2
            # Rewrite `- key: v` as a plain line at the item's own indent so
            # the mapping/scalar logic below handles it uniformly.
            lines[i] = ' ' * item_indent + rest
            val, i = _parse_block(lines, i, item_indent)
            seq.append(val)
        return seq, i
    if not KEY_RE.match(lines[i].strip()):
        return _strip_trailing_comment(lines[i].strip()), i + 1
    mapping = {}
    while True:
        i = _significant(lines, i)
        if i >= len(lines):
            break
        here = _indent(lines[i])
        if here < cur:
            break
        if here > cur:
            raise YamlSubsetError('unexpected indent at line %d: %r' % (i + 1, lines[i]))
        text = lines[i].strip()
        if text.startswith('- '):
            break
        m = KEY_RE.match(text)
        if not m:
            raise YamlSubsetError('not a mapping entry at line %d: %r' % (i + 1, lines[i]))
        key, raw = m.group(1), (m.group(2) or '').strip()
        i += 1
        if raw in ('|', '|-', '|+', '>', '>-'):
            mapping[key], i = _read_block_scalar(lines, i, cur)
        elif raw == '':
            # A block sequence may sit at its parent key's OWN indent
            # (`steps:` / `- uses:` both at 4 here); a nested mapping may not.
            j = _significant(lines, i)
            same_indent_seq = (j < len(lines) and _indent(lines[j]) == cur
                               and lines[j].lstrip().startswith('- '))
            mapping[key], i = _parse_block(lines, i, cur if same_indent_seq else cur + 1)
        else:
            mapping[key] = _strip_trailing_comment(raw)
    return mapping, i


def parse_workflow(path):
    doc, _ = _parse_block(open(path, errors='replace').read().splitlines(), 0, 0)
    if not isinstance(doc, dict) or not isinstance(doc.get('jobs'), dict):
        raise YamlSubsetError('%s has no top-level `jobs:` mapping after parsing '
                              '-- the parser and the file disagree' % path)
    for name, job in doc['jobs'].items():
        if not isinstance(job, dict) or not isinstance(job.get('steps'), list):
            raise YamlSubsetError('job %r has no `steps:` list after parsing' % name)
    return doc


# ----------------------------------------------------------------- checks ---

def invocation_lines(run_body):
    """Lines of a `run:` body that invoke the checker, shell comments removed
    first -- the real step's body is full of `# ...` prose that names it."""
    hits = []
    for raw in run_body.splitlines():
        line = raw.split('#', 1)[0].strip()
        if CHECKER_REL in line or 'check_board_separation.py' in line:
            hits.append(line)
    return hits


def main():
    print('Regression guard: CI must still run check_board_separation.py, at '
          'fetch-depth 0 (rule 21c)')
    problems = []
    doc = parse_workflow(WORKFLOW)

    found = []          # (job_name, step, invocation lines)
    for job_name, job in doc['jobs'].items():
        for step in job['steps']:
            if not isinstance(step, dict):
                continue
            body = step.get('run')
            if not isinstance(body, str):
                continue
            hits = invocation_lines(body)
            if hits:
                found.append((job_name, step, hits))

    print('  parsed %d job(s) from %s' % (len(doc['jobs']), WORKFLOW))
    if not found:
        print('\nFAIL: no step in %s runs `python3 %s`.' % (WORKFLOW, CHECKER_REL))
        print('  The merge-boundary half of rule 21c is then enforced by')
        print('  nothing: testsys/hooks/pre-commit cannot see `git commit')
        print('  --amend`, so a mixed board/code commit lands green. Restore')
        print('  the step, or delete this guard in the same commit that')
        print('  records what replaced it.')
        return 1
    if len(found) > 1:
        problems.append('%d steps invoke the checker (%s) -- one authoritative '
                        'invocation, or a reviewer cannot tell which range is '
                        'actually gated'
                        % (len(found), ', '.join(j for j, _, _ in found)))

    job_name, step, hits = found[0]
    print('  step %r in job %r invokes it: %r'
          % (step.get('name', '<unnamed>'), job_name, hits))

    with_arg = [h for h in hits if re.search(r'check_board_separation\.py\s+\S', h)]
    if not with_arg:
        problems.append('the invocation %r passes NO range argument; the '
                        'checker refuses with exit 2 and gates nothing' % hits)

    if 'if' in step:
        problems.append("the step carries `if: %s` -- a condition can disarm it "
                        "while leaving it visible and the job green"
                        % step['if'])

    job = doc['jobs'][job_name]
    checkouts = [s for s in job['steps']
                 if isinstance(s, dict) and str(s.get('uses', '')).startswith('actions/checkout')]
    if not checkouts:
        problems.append('job %r has no actions/checkout step, so the range the '
                        'checker resolves comes from an unknown tree' % job_name)
    for co in checkouts:
        depth = (co.get('with') or {}).get('fetch-depth')
        if str(depth) != '0':
            problems.append(
                "job %r checks out with fetch-depth=%s (absent means the "
                "action's default, 1). The checker resolves a COMMIT RANGE: at "
                "depth 1 the range endpoints are not in the clone, so the step "
                "cannot check the commits it claims to check."
                % (job_name, depth if depth is not None else 'absent'))

    if not os.path.isfile(CHECKER):
        problems.append('%s does not exist -- the step invokes nothing' % CHECKER)
    else:
        if not os.access(CHECKER, os.X_OK):
            problems.append('%s is not executable (rule 13)' % CHECKER)
        r = subprocess.run([sys.executable, CHECKER], capture_output=True,
                           text=True, timeout=60)
        if r.returncode != 2:
            problems.append(
                'running the checker with no range exited %d, not the 2 its '
                'own contract documents ("exactly one commit range is '
                'required; there is no default") -- it is present but no '
                'longer behaves as the step assumes. stderr: %r'
                % (r.returncode, r.stderr.strip()[:200]))
        else:
            print('  %s: present, executable, refuses a missing range with '
                  'exit 2' % CHECKER_REL)

    if problems:
        print('\nFAIL: %d problem(s):' % len(problems))
        for p in problems:
            print('  - %s' % p)
        return 1
    print('\nPASS: one step (%r, job %r) runs the checker with a range, '
          'unconditionally, at fetch-depth 0; the checker is executable and '
          'refuses an empty range' % (step.get('name', '<unnamed>'), job_name))
    return 0


if __name__ == '__main__':
    sys.exit(main())

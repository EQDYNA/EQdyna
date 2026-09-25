#! /usr/bin/env python3
"""
Regression guard for board row 87 (owner ruling 2026-09-24): test.yml's
`push` trigger runs for master and tags only; `pull_request` stays. A raw
feature-branch push duplicated its PR run (rule 25 sends every branch
through a PR).

Asserted on the parsed `on:` block of .github/workflows/test.yml:
  - push.branches is exactly [master];
  - push.tags is present and non-empty (naming only branches: would stop
    tag pushes from triggering at all);
  - pull_request is still a trigger.
Both ways (rule 14a): the checker first runs on the pre-row-87 trigger text
(push with no branch filter) and must reject it.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
WORKFLOW = os.path.join(ROOT, '.github', 'workflows', 'test.yml')

PRE_ROW_87 = """on:
  push:
    paths-ignore:
      - 'docs/**'
  pull_request:
jobs:
"""


def on_block(text):
    """{trigger: {key: [items]}} for the top-level `on:` block. Indentation
    parse of the shape this workflow uses; no YAML dependency."""
    lines = text.splitlines()
    try:
        start = next(i for i, l in enumerate(lines) if l.rstrip() == 'on:')
    except StopIteration:
        raise ValueError('no top-level `on:` line')
    out, trig, key = {}, None, None
    for l in lines[start + 1:]:
        if not l.strip() or l.lstrip().startswith('#'):
            continue
        ind = len(l) - len(l.lstrip())
        if ind == 0:
            break
        body = l.strip()
        if ind == 2:
            trig = body.rstrip(':'); out[trig] = {}; key = None
        elif ind == 4 and body.endswith(':'):
            key = body.rstrip(':'); out[trig][key] = []
        elif body.startswith('- ') and trig is not None and key is not None:
            out[trig][key].append(body[2:].strip().strip("'\""))
    return out


def problems(text):
    on = on_block(text)
    bad = []
    push = on.get('push')
    if push is None:
        bad.append('no push trigger (master pushes must still run: pr-policy-gate fires only there)')
    else:
        if push.get('branches') != ['master']:
            bad.append('push.branches is %r, want exactly [master]' % push.get('branches'))
        if not push.get('tags'):
            bad.append('push.tags missing/empty -- a branches-only filter stops tag pushes triggering')
    if 'pull_request' not in on:
        bad.append('pull_request trigger removed')
    return bad


def main():
    if not problems(PRE_ROW_87):
        print('FAIL test_ci_push_trigger_filter: negative control (unfiltered push) '
              'was accepted -- the checker cannot go red')
        return 1
    with open(WORKFLOW) as fh:
        bad = problems(fh.read())
    if bad:
        print('FAIL test_ci_push_trigger_filter')
        for b in bad:
            print(' -', b)
        return 1
    print('SUCCESS test_ci_push_trigger_filter: push = master + tags, pull_request '
          'kept; negative control (unfiltered push) rejected')
    return 0


if __name__ == '__main__':
    sys.exit(main())

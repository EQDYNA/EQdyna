#! /usr/bin/env python3
"""testsys/needs_sweep.py <base>..<head> -- item 2, owner 2026-09-30
("fewer full sweeps and fewer releases"): an ADVISORY tool reporting
whether a git range touches any path testsys/change_class.py classifies
PHYSICS, i.e. whether `testsys/run.py e2e` (or `release`) needs to run
again before anyone trusts this range did not change a gated result.

    python3 testsys/needs_sweep.py origin/master..HEAD
    python3 testsys/needs_sweep.py --strict origin/master..HEAD

Prints exactly one line:
  "FULL SWEEP NEEDED: <path>, physics path changed (N physics path(s) of M changed)"
  "fast tier only: no physics path changed (N files)"
Exits 0 either way UNLESS --strict, in which case a needed sweep exits 1 --
for an agent or a CI step that wants a hard gate rather than an advisory
read. Unresolvable input (a bad range spec) fails toward FULL SWEEP NEEDED
(rule 2), same shape as check_release_due.py's own error handling, and
--strict then exits 1 too.

This does NOT consult docs/release_exempt.txt: review item 7 -- an
exemption waives rule 27's release-due-at-once trigger only, never this
tool's sweep verdict. A change exempted as "output-neutral" still needs
the sweep that PROVED it was output-neutral in the first place.
"""
import argparse
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from testsys import change_class  # noqa: E402


def git(*args, cwd=None):
    """cwd defaults to the PROCESS's actual working directory (os.getcwd()),
    never the static location of this file (ROOT): needs_sweep.py is meant
    to be run from inside whatever repo/range it is asked about (its own
    CLI tests invoke it via subprocess against a throwaway tiny_repo), and
    ROOT is only where `testsys` lives for the `import change_class` above,
    not necessarily the repo the caller wants diffed."""
    return subprocess.run(['git', *args], cwd=cwd or os.getcwd(),
                          capture_output=True, text=True, errors='replace',
                          check=True).stdout


def changed_paths(range_spec, cwd=None):
    """Every path touched anywhere across `range_spec` (a two-dot
    `base..head` git range), resolved against `cwd` (default: the caller's
    own working directory -- see git() above)."""
    # --no-renames (same guard as pr_policy.commit_files): with rename
    # detection on, a PHYSICS file moved to an INTERNAL path reports only
    # its NEW name and the range would read as fast-tier-only.
    out = git('diff', '--name-only', '--no-renames', range_spec, cwd=cwd)
    return sorted(set(p for p in out.splitlines() if p.strip()))


def assess(paths):
    """(affects_gate, ships_to_user, output_change, [(path, bucket), ...])
    -- item 6's three flags, computed ONCE here so callers (check_release_
    due.py's own per-commit loop aside, which needs finer per-diff detail)
    can reuse this range-level verdict rather than re-deriving it.

    affects_gate  : classify_path(p) == PHYSICS for some changed p --
                    "needs a sweep" (this tool's own verdict).
    ships_to_user : classify_path(p) == USER_FACING for some changed p.
    output_change : change_class.is_output_change(p) for some changed p --
                    rule 27's narrower "release due at once" criterion
                    (checked WITHOUT a diff here -- a range-level, not a
                    per-commit, view -- so a CODE_DIRS path counts IN
                    unconditionally, same conservative default
                    is_output_change itself documents for diff_text=None).
    """
    classified = [(p, change_class.classify_path(p)) for p in paths]
    affects_gate = any(c == change_class.PHYSICS for _, c in classified)
    ships_to_user = any(c == change_class.USER_FACING for _, c in classified)
    output_change = any(change_class.is_output_change(p) for p in paths)
    return affects_gate, ships_to_user, output_change, classified


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('range', help='a git range, e.g. origin/master..HEAD')
    ap.add_argument('--strict', action='store_true',
                    help='exit 1 (instead of 0) when a full sweep is needed')
    args = ap.parse_args(argv)

    if '..' not in args.range:
        print('needs_sweep: %r is not a two-dot git range (base..head) -- '
              'failing toward FULL SWEEP NEEDED (rule 2)' % args.range)
        return 1 if args.strict else 0
    try:
        paths = changed_paths(args.range)
    except subprocess.CalledProcessError as exc:
        print('needs_sweep: could not resolve %r (%s) -- failing toward '
              'FULL SWEEP NEEDED (rule 2)' % (args.range, exc.stderr.strip()
                                              if exc.stderr else exc))
        return 1 if args.strict else 0

    if not paths:
        print('fast tier only: no physics path changed (0 files)')
        return 0

    affects_gate, ships_to_user, output_change, classified = assess(paths)
    if affects_gate:
        first = next(p for p, c in classified if c == change_class.PHYSICS)
        n_physics = sum(1 for _, c in classified if c == change_class.PHYSICS)
        print('FULL SWEEP NEEDED: %s, physics path changed (%d physics '
              'path(s) of %d changed)' % (first, n_physics, len(paths)))
        return 1 if args.strict else 0

    print('fast tier only: no physics path changed (%d files)' % len(paths))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

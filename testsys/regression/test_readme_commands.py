#! /usr/bin/env python3
"""
Regression guard: commands README tells a user to run must exist (rules 2, 11).

THREE DOCUMENTED COMMANDS HAVE BEEN WRONG IN ONE DAY:

  1. `python -m eqdyna.standalone` (README's GPU line) -- the module was
     deleted with the standalone package. Raised ModuleNotFoundError.
  2. `python3 testsys/run.py e2e-full` -- real, but it REFUSES to launch
     without EQDYNA_FULL_LAUNCH=yes-hours, which README did not mention, so a
     user following it verbatim got `FAIL e2e-full (exit 1)`.
  3. PROJECT_RULES rule 16 told the reader to reproduce CI's entry point and
     then quoted `run.py all`, while CI runs `unit regression e2e-ci`.

Each was a doc that had drifted from a tree that moved underneath it. A README
command that does not work is worse than an undocumented one: the reader
trusts it and concludes the software is broken.

WHAT THIS PINS, and deliberately not more. It checks that every command README
presents as runnable RESOLVES -- the script exists, the tier name is real, the
module is importable, the referenced file is present. It does NOT execute the
expensive ones; running the sweep from a regression guard would make the cheap
tier cost 20 minutes (rule 9).

So this catches the class "documented thing does not exist" and not the class
"documented thing exists but misbehaves". The second needs a real run, which
is what the e2e tier is for.

Cheap (rule 9): filesystem and import checks only, under 2 s.
Exits non-zero on any failure.
"""
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
README = os.path.join(ROOT, 'README.md')


def readme():
    return open(README, errors='replace').read()


def check_referenced_scripts_exist():
    """Every repo-relative path README shows in a command must be present."""
    text = readme()
    # paths that look like a script or module invocation in a command position
    paths = set(re.findall(r'(?:python3?\s+|bash\s+|\./)([\w./-]+\.(?:py|sh))', text))
    paths |= set(re.findall(r'`(testsys/[\w./-]+\.py)`', text))
    # Generated INSIDE a case directory by create.newcase/case.setup, so they
    # are not repo-relative paths and must not be looked for at the root.
    paths -= {'run.sh', 'case.setup', 'case.submit', 'clean.py'}
    missing = [p for p in sorted(paths)
               if not os.path.exists(os.path.join(ROOT, p))]
    if missing:
        raise AssertionError('README names script(s) that do not exist: %r' % missing)
    print('  PASS  %d referenced script path(s) all exist' % len(paths))


def check_tier_names_are_real():
    """`testsys/run.py <tier>` -- every tier README names must be known to
    run.py, or the user gets 'unknown tier'."""
    text = readme()
    named = set(re.findall(r'run\.py\s+([a-z0-9\-]+)', text))
    named.discard('py')
    sys.path.insert(0, os.path.join(ROOT, 'testsys'))
    import importlib.util as u
    spec = u.spec_from_file_location('runpy_mod', os.path.join(ROOT, 'testsys', 'run.py'))
    mod = u.module_from_spec(spec)
    spec.loader.exec_module(mod)
    known = set(mod.TIERS) | set(mod.OPTIONAL_TIERS) | {'all'}
    unknown = sorted(named - known)
    if unknown:
        raise AssertionError(
            'README names tier(s) run.py does not know: %r (known: %s)'
            % (unknown, ', '.join(sorted(known))))
    print('  PASS  %d tier name(s) named by README all exist: %s'
          % (len(named), ', '.join(sorted(named))))


def check_python_module_invocations_import():
    """`python3 -m <mod>` must be importable with src/python on the path.

    This is the check that would have caught `python -m eqdyna.standalone`.
    """
    text = readme()
    mods = set(re.findall(r'python3?\s+-m\s+([\w.]+)', text))
    if not mods:
        print('  PASS  README invokes no python -m module (nothing to check)')
        return
    env = dict(os.environ)
    env['PYTHONPATH'] = os.pathsep.join(
        [os.path.join(ROOT, 'src', 'python'), env.get('PYTHONPATH', '')])
    bad = []
    for m in sorted(mods):
        r = subprocess.run([sys.executable, '-c', 'import %s' % m],
                           env=env, capture_output=True, text=True)
        if r.returncode != 0:
            bad.append((m, r.stderr.strip().splitlines()[-1] if r.stderr else '?'))
    if bad:
        raise AssertionError(
            'README invokes python -m on module(s) that do not import: %s'
            % '; '.join('%s (%s)' % b for b in bad))
    print('  PASS  python -m module(s) importable: %s' % ', '.join(sorted(mods)))


def check_optin_env_vars_are_documented():
    """A tier that refuses without an env var must have that var in README
    beside it, or the documented command simply fails."""
    text = readme()
    problems = []
    for script, var in (('testsys/e2e/run_e2e_full.py', 'EQDYNA_FULL_LAUNCH'),):
        p = os.path.join(ROOT, script)
        if not os.path.exists(p):
            continue
        if var not in open(p, errors='replace').read():
            continue                       # the script no longer gates on it
        if var not in text:
            problems.append(
                '%s refuses to run without %s, and README never mentions it'
                % (script, var))
    if problems:
        raise AssertionError('; '.join(problems))
    print('  PASS  opt-in environment variables are documented in README')


def check_compsets_named_by_readme_exist():
    text = readme()
    names = set(re.findall(r'case_input/(test\.[\w.]+)', text))
    names |= set(re.findall(r'^\* (test\.[\w.]+)', text, re.M))
    missing = [n for n in sorted(names)
               if not os.path.isdir(os.path.join(ROOT, 'case_input', n))]
    if missing:
        raise AssertionError('README names compset(s) that do not exist: %r' % missing)
    print('  PASS  %d compset(s) named by README all exist' % len(names))


def main():
    print('Regression guard: README commands must resolve')
    checks = [check_referenced_scripts_exist,
              check_tier_names_are_real,
              check_python_module_invocations_import,
              check_optin_env_vars_are_documented,
              check_compsets_named_by_readme_exist]
    failures = []
    for c in checks:
        try:
            c()
        except AssertionError as e:
            failures.append(str(e))
            print('  FAIL  %s: %s' % (c.__name__, e))
    if failures:
        print('\nFAIL test_readme_commands (%d check(s))' % len(failures))
        return 1
    print('\nSUCCESS test_readme_commands')
    return 0


if __name__ == '__main__':
    sys.exit(main())

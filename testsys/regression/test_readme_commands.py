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


def check_referenced_scripts_exist(text=None, label='README.md'):
    """Every repo-relative path README shows in a command must be present.

    Parameterized (text/label) so a second guard covering docs/user/**/*.md
    can call this same resolution logic against a different document instead
    of re-implementing it; `text=None` (README.md's own use) is unchanged."""
    text = readme() if text is None else text
    # paths that look like a script or module invocation in a command position
    paths = set(re.findall(r'(?:python3?\s+|bash\s+|\./)([\w./-]+\.(?:py|sh))', text))
    paths |= set(re.findall(r'`(testsys/[\w./-]+\.py)`', text))
    # Generated INSIDE a case directory by create.newcase/case.setup, so they
    # are not repo-relative paths and must not be looked for at the root.
    paths -= {'run.sh', 'case.setup', 'case.submit', 'clean.py'}
    missing = [p for p in sorted(paths)
               if not os.path.exists(os.path.join(ROOT, p))]
    if missing:
        raise AssertionError('%s names script(s) that do not exist: %r' % (label, missing))
    print('  PASS  %s: %d referenced script path(s) all exist' % (label, len(paths)))


def check_tier_names_are_real(text=None, label='README.md'):
    """`testsys/run.py <tier>` -- every tier a document names must be known
    to run.py, or the user gets 'unknown tier'."""
    text = readme() if text is None else text
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
            '%s names tier(s) run.py does not know: %r (known: %s)'
            % (label, unknown, ', '.join(sorted(known))))
    print('  PASS  %s: %d tier name(s) all exist: %s'
          % (label, len(named), ', '.join(sorted(named))))


def check_python_module_invocations_import(text=None, label='README.md'):
    """`python3 -m <mod>` must be importable with src/python on the path.

    This is the check that would have caught `python -m eqdyna.standalone`.
    """
    text = readme() if text is None else text
    mods = set(re.findall(r'python3?\s+-m\s+([\w.]+)', text))
    # `pip` is the installer itself, and a module the SAME document installs
    # with `pip install` (e.g. virtualenv) cannot import before its own
    # install step runs -- test_readme_executes proves those by running them.
    # Every other module must import here.
    installed = set()
    for line in re.findall(r'pip3?\s+install\s+([^\n#]+)', text):
        for tok in line.split():
            if not tok.startswith('-'):
                installed.add(re.split(r'[\[<>=]', tok.strip('"\''))[0].lower())
    mods = {m for m in mods if m != 'pip' and m.split('.')[0].lower() not in installed}
    if not mods:
        print('  PASS  %s: invokes no python -m module (nothing to check)' % label)
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
            '%s invokes python -m on module(s) that do not import: %s'
            % (label, '; '.join('%s (%s)' % b for b in bad)))
    print('  PASS  %s: python -m module(s) importable: %s' % (label, ', '.join(sorted(mods))))


def check_optin_env_vars_are_documented(text=None, label='README.md'):
    """A tier that refuses without an env var must have that var documented
    beside it, or the documented command simply fails.

    Only enforced on a document that actually invokes the gated script (its
    basename) or names its opt-in tier ('e2e-full') -- reused across many
    docs/user pages, most of which never mention e2e-full at all, and a page
    that never tells a user to run the script has nothing to fail here."""
    text = readme() if text is None else text
    problems = []
    for script, var in (('testsys/e2e/run_e2e_full.py', 'EQDYNA_FULL_LAUNCH'),):
        p = os.path.join(ROOT, script)
        if not os.path.exists(p):
            continue
        if var not in open(p, errors='replace').read():
            continue                       # the script no longer gates on it
        if os.path.basename(script) not in text and 'e2e-full' not in text:
            continue                       # this doc never invokes it
        if var not in text:
            problems.append(
                '%s refuses to run without %s, and %s never mentions it'
                % (script, var, label))
    if problems:
        raise AssertionError('; '.join(problems))
    print('  PASS  %s: opt-in environment variables are documented' % label)


def check_compsets_named_by_readme_exist(text=None, label='README.md'):
    text = readme() if text is None else text
    names = set(re.findall(r'case_input/(test\.[\w.]+)', text))
    names |= set(re.findall(r'^\* (test\.[\w.]+)', text, re.M))
    missing = [n for n in sorted(names)
               if not os.path.isdir(os.path.join(ROOT, 'case_input', n))]
    if missing:
        raise AssertionError('%s names compset(s) that do not exist: %r' % (label, missing))
    print('  PASS  %s: %d compset(s) named all exist' % (label, len(names)))


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

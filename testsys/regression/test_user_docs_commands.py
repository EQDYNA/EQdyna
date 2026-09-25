#! /usr/bin/env python3
"""
Regression guard: every command shown on the user docs site (docs/user/**/*.md)
must resolve -- the same contract test_readme_commands.py already holds
README.md to (rules 2, 11), extended to the new doc tree (board row 126).

Reuses test_readme_commands.py's own check functions (parameterized there by
text/label) rather than re-implementing script/tier/module/compset
resolution a second time. This file adds no new resolution logic of its own;
it only supplies a different (text, label) per doc page under docs/user/.

Cheap (rule 9): filesystem, import and one importlib load per doc page,
under a couple of seconds for the whole site.
Exits non-zero on any failure.
"""
import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import test_readme_commands as trc  # noqa: E402  (path insert must come first)

ROOT = trc.ROOT


def user_doc_files():
    return sorted(
        os.path.relpath(p, ROOT)
        for p in glob.glob(os.path.join(ROOT, 'docs', 'user', '**', '*.md'), recursive=True)
    )


def main():
    print('Regression guard: docs/user/**/*.md commands must resolve')
    files = user_doc_files()
    if not files:
        print('FAIL test_user_docs_commands: no docs/user/**/*.md files found')
        return 1
    checks = [
        trc.check_referenced_scripts_exist,
        trc.check_tier_names_are_real,
        trc.check_python_module_invocations_import,
        trc.check_optin_env_vars_are_documented,
        trc.check_compsets_named_by_readme_exist,
    ]
    failures = []
    for rel in files:
        text = open(os.path.join(ROOT, rel), errors='replace').read()
        for fn in checks:
            try:
                fn(text, rel)
            except AssertionError as e:
                failures.append(str(e))
                print('  FAIL  %s (%s): %s' % (fn.__name__, rel, e))
    if failures:
        print('\nFAIL test_user_docs_commands (%d check(s), %d doc file(s))'
              % (len(failures), len(files)))
        return 1
    print('\nSUCCESS test_user_docs_commands: %d doc file(s) checked' % len(files))
    return 0


if __name__ == '__main__':
    sys.exit(main())

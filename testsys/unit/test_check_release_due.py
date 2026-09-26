"""Unit tests for testsys/regression/check_release_due.py (rule 27, item 133).

The classifier must go both ways (rule 14a): a real code change counts, and
a comment-only change does not. Anything it cannot classify counts IN.
"""
import importlib.util
import os
import subprocess
import sys

from conftest import REPO_ROOT

_spec = importlib.util.spec_from_file_location(
    'check_release_due',
    os.path.join(REPO_ROOT, 'testsys', 'regression', 'check_release_due.py'))
crd = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(crd)


def _diff(minus, plus):
    return ('diff --git a/x b/x\n--- a/x\n+++ b/x\n@@ -1 +1 @@\n'
            + ''.join('-%s\n' % l for l in minus)
            + ''.join('+%s\n' % l for l in plus))


def test_fortran_code_change_counts():
    assert crd.classify('src/fortran/fric.f90', _diff(['x = 1'], ['x = 2']))


def test_fortran_comment_only_change_does_not_count():
    assert crd.classify('src/fortran/fric.f90',
                        _diff(['! old note'], ['! new note', ''])) is None


def test_python_comment_only_change_does_not_count():
    assert crd.classify('src/python/eqdyna/fric.py',
                        _diff(['# old'], ['    # new'])) is None


def test_python_code_change_counts():
    assert crd.classify('src/python/eqdyna/fric.py', _diff(['a = 1'], ['a = 3']))


def test_docstring_change_counts_in_conservatively():
    assert crd.classify('src/python/eqdyna/fric.py',
                        _diff(['    """Old words."""'], ['    """New words."""']))


def test_output_writer_counts_even_if_comment_only():
    assert crd.classify('src/fortran/library_output.f90', _diff(['! a'], ['! b']))


def test_reference_change_counts():
    assert crd.classify('test.reference.results/test.tpv8/frt.canonical.txt', '')


def test_unclassifiable_code_diff_counts_in():
    assert crd.classify('src/fortran/fric.f90', '')


def test_docs_change_does_not_count():
    assert crd.classify('docs/user/benchmarks.md', _diff(['a'], ['b'])) is None


def test_script_prints_one_line_and_exits_zero():
    r = subprocess.run([sys.executable, os.path.join(
        REPO_ROOT, 'testsys', 'regression', 'check_release_due.py')],
        capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    lines = [l for l in r.stdout.splitlines() if l.strip()]
    assert len(lines) == 1, r.stdout
    assert lines[0].startswith(('RELEASE DUE:', 'release not due:')), lines[0]

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


def _run_main(monkeypatch, capsys, days, subjects, files, tmp_path=None,
              exempt_lines=()):
    """`files` is the list of paths changed by ONE commit since the tag
    (sha c0ffee1...); `exempt_lines` are written to a temp exemption file."""
    import datetime
    tag_date = (datetime.datetime.now(datetime.timezone.utc)
                - datetime.timedelta(days=days)).isoformat()
    sha = 'c0ffee1' + '0' * 33

    def fake_git(*a):
        if a[0] == 'describe':
            return 'vX\n'
        if a[0] == 'log' and '--format=%cI' in a:
            return tag_date + '\n'
        if a[0] == 'log' and '--format=%H' in a:
            return (sha + '\n') if files else ''
        if a[0] == 'log':
            return subjects
        if a[0] == 'diff-tree':
            return '\n'.join(files) + '\n'
        return _diff(['x = 1'], ['x = 2'])
    monkeypatch.setattr(crd, 'git', fake_git)
    path = os.path.join(str(tmp_path or '/nonexistent-dir'), 'exempt.txt')
    if tmp_path is not None:
        open(path, 'w').write(''.join(l + '\n' for l in exempt_lines))
    monkeypatch.setattr(crd, 'EXEMPT_PATH', path)
    assert crd.main() == 0
    return capsys.readouterr().out.strip()


CODE = ['src/fortran/fric.f90']
FIVE = ''.join('fix (#%d)\n' % n for n in range(5))
FOUR = ''.join('fix (#%d)\n' % n for n in range(4))


def test_physics_change_is_due_at_once(monkeypatch, capsys):
    """Owner's wording: 'as soon as a physics or output change lands' --
    0 days, 1 PR is already due."""
    out = _run_main(monkeypatch, capsys, 0, 'x (#1)\n', CODE)
    assert out.startswith('RELEASE DUE: physics/output change'), out


def test_threshold_is_due_without_a_physics_change(monkeypatch, capsys):
    """'or at the latest after a week or ~5 PRs' (owner lowered 10 -> 5,
    2026-09-28): PRs that are not physics still make a release due once 7
    days or 5 PRs have passed, and 4 PRs inside a week does not."""
    assert _run_main(monkeypatch, capsys, 7, 'x (#1)\n',
                     ['docs/a.md']).startswith('RELEASE DUE: the days threshold')
    assert _run_main(monkeypatch, capsys, 0, FOUR,
                     ['docs/a.md']).startswith('release not due')
    assert _run_main(monkeypatch, capsys, 0, FIVE,
                     ['docs/a.md']).startswith('RELEASE DUE: the PRs threshold')
    assert _run_main(monkeypatch, capsys, 6, 'x (#1)\n',
                     ['docs/a.md']).startswith('release not due')


def test_docs_only_stretch_with_no_pr_never_forces(monkeypatch, capsys):
    assert _run_main(monkeypatch, capsys, 30, 'board: x\n',
                     ['docs/a.md']).startswith('release not due')


def test_exempt_commit_does_not_trigger_but_counts_as_pr(monkeypatch, capsys, tmp_path):
    ev = ['c0ffee1 bit-identical: run.py all 23/23']
    out = _run_main(monkeypatch, capsys, 0, 'x (#1)\n', CODE, tmp_path, ev)
    assert out.startswith('release not due') and '1 exempt' in out, out
    out = _run_main(monkeypatch, capsys, 7, 'x (#1)\n', CODE, tmp_path, ev)
    assert out.startswith('RELEASE DUE: the days threshold'), out


def test_exemption_without_evidence_exempts_nothing(monkeypatch, capsys, tmp_path):
    for bad in (['c0ffee1'], ['c0ffee1   '], ['zzzzzzz evidence'], ['c0ff evidence']):
        out = _run_main(monkeypatch, capsys, 0, 'x (#1)\n', CODE, tmp_path, bad)
        assert out.startswith('RELEASE DUE: physics/output change'), (bad, out)


def test_merge_commit_subjects_count_as_prs(monkeypatch, capsys):
    subs = ''.join('Merge pull request #%d from x/y\n' % n for n in range(5))
    assert _run_main(monkeypatch, capsys, 0, subs,
                     ['docs/a.md']).startswith('RELEASE DUE: the PRs threshold')


def test_git_error_prints_one_due_line_and_exits_zero(monkeypatch, capsys):
    def boom(*a):
        if a[0] == 'describe':
            return 'vX\n'
        raise subprocess.CalledProcessError(128, ['git'] + list(a),
                                            stderr='fatal: bad object')
    monkeypatch.setattr(crd, 'git', boom)
    assert crd.main() == 0
    out = capsys.readouterr().out.strip().splitlines()
    assert len(out) == 1 and out[0].startswith('RELEASE DUE: checker error'), out


def test_non_utf8_fortran_comment_does_not_crash(monkeypatch, capsys, tmp_path):
    """A Latin-1 comment with no NUL byte is not 'Binary' to git but is not
    UTF-8: decoding it strictly used to raise (PR #49 audit, measured). A
    real throwaway repo, so the real git() decode path runs."""
    def g(*a):
        subprocess.run(['git', *a], cwd=tmp_path, check=True,
                       capture_output=True)
    g('init', '-q'); g('config', 'user.email', 't@t'); g('config', 'user.name', 't')
    (tmp_path / 'README').write_text('x\n'); g('add', '.'); g('commit', '-qm', 'root')
    g('tag', 'v0')
    f = tmp_path / 'src' / 'fortran'; f.mkdir(parents=True)
    (f / 'fric.f90').write_bytes(b'! caf\xe9 note\nx = 1\n')
    g('add', '.'); g('commit', '-qm', 'latin1 fortran (#1)')
    monkeypatch.setattr(crd, 'ROOT', str(tmp_path))
    monkeypatch.setattr(crd, 'EXEMPT_PATH', str(tmp_path / 'none.txt'))
    assert crd.main() == 0
    out = capsys.readouterr().out.strip().splitlines()
    assert len(out) == 1, out
    assert out[0].startswith('RELEASE DUE: physics/output change'), out

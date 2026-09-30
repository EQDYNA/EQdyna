"""Unit tests for testsys/needs_sweep.py (item 2, owner 2026-09-30).

Both ways (rule 14a): a physics path in the range prints FULL SWEEP NEEDED
and (with --strict) exits 1; a range with none prints the fast-tier line
and always exits 0.
"""
import os
import sys
import subprocess

import pytest

from conftest import REPO_ROOT
from testsys import change_class, needs_sweep

NEEDS_SWEEP = os.path.join(REPO_ROOT, 'testsys', 'needs_sweep.py')


# --------------------------------------------------------------------------
# assess(): the range-level classification, fixture-driven (Victor review's
# named cases): docs no; testNameList.py yes; scripts/lib.py yes;
# perf/ledger.py yes; the classifier file itself yes.
# --------------------------------------------------------------------------
@pytest.mark.parametrize('paths,want_affects_gate', [
    (['docs/notes/foo.md'], False),
    (['pathway_forward.md'], False),
    (['testsys/regression/test_something.py'], False),
    (['testNameList.py'], True),
    (['scripts/lib.py'], True),
    (['testsys/perf/ledger.py'], True),
    (['testsys/change_class.py'], True),
    (['src/fortran/driver.f90'], True),
    (['case_input/test.tpv8/user_defined_params.py'], True),
])
def test_assess_affects_gate_per_fixture(paths, want_affects_gate):
    affects_gate, _, _, _ = needs_sweep.assess(paths)
    assert affects_gate is want_affects_gate, (paths, affects_gate)


def test_assess_ships_to_user_for_docs_user_and_readme():
    affects_gate, ships_to_user, _, _ = needs_sweep.assess(['docs/user/getting-started.md'])
    assert not affects_gate
    assert ships_to_user


def test_assess_mixed_range_any_physics_path_wins():
    affects_gate, ships_to_user, _, classified = needs_sweep.assess(
        ['docs/notes/foo.md', 'testsys/matrix.py', 'README.md'])
    assert affects_gate
    assert ships_to_user
    assert dict(classified)['testsys/matrix.py'] == change_class.PHYSICS


# --------------------------------------------------------------------------
# The CLI, against a real throwaway git repo.
# --------------------------------------------------------------------------
def _git(tmp_path, *args):
    subprocess.run(['git', *args], cwd=tmp_path, check=True, capture_output=True)


@pytest.fixture
def tiny_repo(tmp_path):
    _git(tmp_path, 'init', '-q')
    _git(tmp_path, 'config', 'user.email', 't@t')
    _git(tmp_path, 'config', 'user.name', 't')
    (tmp_path / 'README.md').write_text('x\n')
    _git(tmp_path, 'add', '.')
    _git(tmp_path, 'commit', '-qm', 'root')
    _git(tmp_path, 'branch', 'base')
    return tmp_path


def test_cli_docs_only_range_is_fast_tier(tiny_repo):
    (tiny_repo / 'docs').mkdir()
    (tiny_repo / 'docs' / 'note.md').write_text('note\n')
    _git(tiny_repo, 'add', '.')
    _git(tiny_repo, 'commit', '-qm', 'docs')
    r = subprocess.run([sys.executable, NEEDS_SWEEP, 'base..HEAD'],
                       cwd=tiny_repo, capture_output=True, text=True)
    assert r.returncode == 0
    assert r.stdout.strip().startswith('fast tier only'), r.stdout


def test_cli_physics_path_needs_sweep_and_strict_exits_1(tiny_repo):
    (tiny_repo / 'testNameList.py').write_text('nameList = []\n')
    _git(tiny_repo, 'add', '.')
    _git(tiny_repo, 'commit', '-qm', 'add a case')
    r = subprocess.run([sys.executable, NEEDS_SWEEP, 'base..HEAD'],
                       cwd=tiny_repo, capture_output=True, text=True)
    assert r.returncode == 0  # advisory by default
    assert r.stdout.strip().startswith('FULL SWEEP NEEDED'), r.stdout
    assert 'testNameList.py' in r.stdout

    r = subprocess.run([sys.executable, NEEDS_SWEEP, '--strict', 'base..HEAD'],
                       cwd=tiny_repo, capture_output=True, text=True)
    assert r.returncode == 1, r.stdout
    assert r.stdout.strip().startswith('FULL SWEEP NEEDED'), r.stdout


def test_cli_bad_range_fails_toward_full_sweep_needed(tiny_repo):
    r = subprocess.run([sys.executable, NEEDS_SWEEP, 'no-such-ref..HEAD'],
                       cwd=tiny_repo, capture_output=True, text=True)
    assert r.returncode == 0
    assert 'FULL SWEEP NEEDED' in r.stdout

    r = subprocess.run([sys.executable, NEEDS_SWEEP, '--strict', 'no-such-ref..HEAD'],
                       cwd=tiny_repo, capture_output=True, text=True)
    assert r.returncode == 1


def test_cli_rename_out_of_physics_path_still_needs_sweep(tiny_repo):
    """PR #62 audit: a src/ file renamed to docs/ must report the OLD path
    (--no-renames); with rename detection it read as fast-tier-only."""
    (tiny_repo / 'src').mkdir()
    (tiny_repo / 'src' / 'a.f90').write_text('x\n')
    _git(tiny_repo, 'add', '.')
    _git(tiny_repo, 'commit', '-qm', 'add src')
    _git(tiny_repo, 'branch', '-f', 'base')
    (tiny_repo / 'docs').mkdir()
    _git(tiny_repo, 'mv', 'src/a.f90', 'docs/a.f90')
    _git(tiny_repo, 'commit', '-qm', 'move out')
    r = subprocess.run([sys.executable, NEEDS_SWEEP, 'base..HEAD'],
                       cwd=tiny_repo, capture_output=True, text=True)
    assert r.stdout.strip().startswith('FULL SWEEP NEEDED'), r.stdout
    assert 'src/a.f90' in r.stdout

"""Board row 91d (owner, 2026-09-26): the perf tools share ONE default
--exclude-cpus set, run_mpi_scaling.DEFAULT_EXCLUDE_CPUS. This checks the
default the tools actually expose, not their source text, by capturing
the parser each tool builds when its main() reaches parse_args.

EQDYNAROOT is pinned to THIS checkout for the duration: run_jaxmpi_ab's
root-must-be-this-checkout gate runs before its parser exists, and a stale
EQDYNAROOT (normal in a worktree) would otherwise exit early and read as a
wrong default instead of the real cause."""
import argparse
import importlib.util
import os
import sys

import pytest

from conftest import REPO_ROOT

PERF = os.path.join(REPO_ROOT, 'testsys', 'perf')
TOOLS = ('run_mpi_scaling', 'run_jaxmpi_ab', 'run_setup_probe')


def _load(name, monkeypatch):
    monkeypatch.setenv('EQDYNAROOT', REPO_ROOT)
    sys.path.insert(0, PERF)
    spec = importlib.util.spec_from_file_location(name, os.path.join(PERF, name + '.py'))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _exclude_default(mod, monkeypatch):
    seen = {}
    real = argparse.ArgumentParser.parse_args

    def spy(self, *a, **k):
        for act in self._actions:
            if '--exclude-cpus' in act.option_strings:
                seen['default'] = act.default
        raise SystemExit(0)
    monkeypatch.setattr(argparse.ArgumentParser, 'parse_args', spy)
    monkeypatch.setattr(sys, 'argv', [mod.__name__])
    with pytest.raises(SystemExit):
        mod.main()
    monkeypatch.setattr(argparse.ArgumentParser, 'parse_args', real)
    assert 'default' in seen, ('%s exited before building its parser; '
                               'the --exclude-cpus default was never read' % mod.__name__)
    return seen['default']


@pytest.mark.parametrize('tool', TOOLS)
def test_every_perf_tool_uses_the_shared_default(tool, monkeypatch):
    shared = _load('run_mpi_scaling', monkeypatch).DEFAULT_EXCLUDE_CPUS
    assert shared == '0,1,16,17,18,19'
    assert _exclude_default(_load(tool, monkeypatch), monkeypatch) == shared

"""Unit tests for create.newcase's case_copy() helper.

This is the function-level counterpart to
testsys/regression/test_create_newcase.py, which guards the same behaviour
at the subprocess/CLI level (the IsADirectoryError incident,
PROJECT_RULES.md rule 10). Here we call case_copy() directly against
throwaway directories: no subprocess, no case_input/, sub-millisecond.
"""
import importlib.machinery
import importlib.util
import os

from conftest import REPO_ROOT


def _load_case_copy():
    # create.newcase has no .py suffix, so spec_from_file_location can't
    # infer a loader on its own -- give it one explicitly.
    path = os.path.join(REPO_ROOT, 'scripts', 'create.newcase')
    loader = importlib.machinery.SourceFileLoader('create_newcase_under_test', path)
    spec = importlib.util.spec_from_loader(loader.name, loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    return mod.case_copy


case_copy = _load_case_copy()


def test_case_copy_copies_plain_files(tmp_path):
    src = tmp_path / 'src'
    dst = tmp_path / 'dst'
    src.mkdir()
    dst.mkdir()
    (src / 'a.txt').write_text('hello')
    (src / 'b.py').write_text('print(1)')

    case_copy(str(src), str(dst))

    assert (dst / 'a.txt').read_text() == 'hello'
    assert (dst / 'b.py').read_text() == 'print(1)'


def test_case_copy_skips_subdirectories(tmp_path):
    # This is the exact shape of the incident being guarded: scripts/ once
    # contained a subdirectory (fractal_stress_diamond_square/) that broke
    # every case creation with an IsADirectoryError.
    src = tmp_path / 'src'
    dst = tmp_path / 'dst'
    src.mkdir()
    dst.mkdir()
    (src / 'keep.py').write_text('ok')
    (src / 'a_tool_collection').mkdir()
    (src / 'a_tool_collection' / 'nested.py').write_text('should not be copied')

    case_copy(str(src), str(dst))  # must not raise IsADirectoryError

    copied = os.listdir(dst)
    assert copied == ['keep.py']


def test_case_copy_on_empty_source_leaves_destination_empty(tmp_path):
    src = tmp_path / 'src'
    dst = tmp_path / 'dst'
    src.mkdir()
    dst.mkdir()

    case_copy(str(src), str(dst))

    assert os.listdir(dst) == []

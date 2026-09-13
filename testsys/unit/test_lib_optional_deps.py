#! /usr/bin/env python3
"""
Regression guard (rule 10, CI failure 2026-09-13): scripts/lib.py must
import WITHOUT optional dependencies (imageio) present — a top-level
import broke every case-runtime script on CI where imageio isn't
installed. Simulates absence via a blocking meta-path finder.
"""
import importlib, importlib.util, sys, os

class _Block:
    blocked = {'imageio'}
    def find_spec(self, name, path=None, target=None):
        if name.split('.')[0] in self.blocked:
            raise ImportError(f'{name} blocked by test (optional dep)')

def test_lib_imports_without_imageio():
    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    blocker = _Block()
    saved = {k: v for k, v in sys.modules.items() if k.split('.')[0] in blocker.blocked}
    for k in saved:
        del sys.modules[k]
    sys.meta_path.insert(0, blocker)
    try:
        spec = importlib.util.spec_from_file_location(
            'lib_isolated', os.path.join(root, 'scripts', 'lib.py'))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)   # must not raise
        assert hasattr(mod, 'loadFrtData') and hasattr(mod, 'generate_gif')
    finally:
        sys.meta_path.remove(blocker)
        sys.modules.update(saved)

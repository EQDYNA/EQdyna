#! /usr/bin/env python3
"""Regression guard: run_e2e._perf_meta must build a snapshot without raising.

Incident (2026-09-23, combo branch a750698): re-applying the profile
collection hooks by hand dropped `_perf_meta`'s local `import ledger`, so
every sweep printed "WARNING: perf-ledger capture failed (NameError: name
'ledger' is not defined)" and appended ZERO ledger rows. `_capture_perf`
deliberately degrades every error to a WARNING (the ledger is never part of a
verdict), so nothing went red: a whole sweep's timings vanished silently.

This calls the real `_perf_meta` exactly the way `_capture_perf` does (same
sys.path insert) on one fake passing cell and asserts a snapshot with that
cell comes back. Mutation: delete the import -> NameError -> RED.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'testsys', 'e2e'))
sys.path.insert(0, os.path.join(ROOT, 'testsys', 'perf'))   # as _capture_perf does

import run_e2e  # noqa: E402


def main():
    results = [('test.tpv8', 'fortran', True, 12.5,
                ['max|diff|=3.051760e-11 bound=1.0e-08'])]
    try:
        meta = run_e2e._perf_meta(results, 'guard', 'cpu', 4, 'deadbee', False)
    except Exception as exc:                        # noqa: BLE001
        print('FAIL test_perf_meta_imports: _perf_meta raised %s: %s -- the '
              'perf ledger would silently get zero rows for every sweep'
              % (type(exc).__name__, exc))
        return 1
    cells = meta.get('cells') if isinstance(meta, dict) else None
    if not cells or len(cells) != 1:
        print('FAIL test_perf_meta_imports: expected 1 cell in the snapshot, '
              'got %r' % (cells,))
        return 1
    print('SUCCESS test_perf_meta_imports (1 fake passing cell -> 1 snapshot cell)')
    return 0


if __name__ == '__main__':
    sys.exit(main())

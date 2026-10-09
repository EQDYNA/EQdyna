#! /usr/bin/env python3
"""
Regression guard (board, 2026-10-09, defect B): TPV12/13 fault-border nodes
must be rupture-eligible, not pinned unbreakable.

TPV12_13_Description_v6.pdf Part 2, "Fault Geometry", p.3: "A node which
lies exactly on the border of the 30000 m x 15000 m rectangle is considered
to be inside the rectangle, and so should be permitted to rupture." Before
this fix, case_input/test.tpv12 and test.tpv13's user_defined_params.py set
mu_s = 1000 (unbreakable) on every along-strike-edge (|x| == fxmax) and
down-dip-edge (z == fzmin) node -- at the gate dx (500 m) that pinned 100%
of the never-ruptured nodes (measured: 182/182 for tpv12, 292/292 for tpv13
unruptured-after-fix nodes include the border set; pre-fix ALL unruptured
nodes were border nodes). This test does not run the solver (rule 9): it
only imports each case's user_defined_params.py (the same module
scripts/case.setup imports) and inspects the mu_s (on_fault_vars[...,1])
array it builds, directly.

Both ways (rule 14a): this FAILS on the pre-fix code (mu_s==1000 on every
border node) and PASSES on the fix (border nodes carry their ordinary
nucleation-zone-weighted mu_s, same as any other node outside the
nucleation patch).
"""
import importlib
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
UNBREAKABLE = 1000.0


def _load_par(case):
    case_dir = os.path.join(ROOT, 'case_input', case)
    saved_path = list(sys.path)
    saved_modules = {k: v for k, v in sys.modules.items()
                      if k in ('user_defined_params', 'defaultParameters', 'lib')}
    for m in ('user_defined_params', 'defaultParameters', 'lib'):
        sys.modules.pop(m, None)
    sys.path.insert(0, case_dir)
    sys.path.insert(0, os.path.join(ROOT, 'scripts'))
    try:
        mod = importlib.import_module('user_defined_params')
        return mod.par
    finally:
        sys.path[:] = saved_path
        for m in ('user_defined_params', 'defaultParameters', 'lib'):
            sys.modules.pop(m, None)
        sys.modules.update(saved_modules)


def _check(case):
    par = _load_par(case)
    mu_s = par.on_fault_vars[:, :, 1]
    fx = par.fx
    fz = par.fz
    fxmax = par.fxmax
    fzmin = fz.min()
    x_border = np.isclose(np.abs(fx), fxmax, atol=0.01)
    z_border = np.isclose(fz, fzmin, atol=0.01)
    border_mask = x_border[np.newaxis, :] | z_border[:, np.newaxis]
    n_border = int(border_mask.sum())
    n_pinned = int((mu_s[border_mask] >= UNBREAKABLE).sum())
    return n_border, n_pinned


def main():
    fails = []
    for case in ('test.tpv12', 'test.tpv13'):
        n_border, n_pinned = _check(case)
        ok = n_border > 0 and n_pinned == 0
        print('  %-12s border nodes=%4d  pinned-unbreakable=%4d  %s'
              % (case, n_border, n_pinned, 'ok' if ok else 'FAIL'))
        if not ok:
            fails.append('%s: %d of %d border nodes still pinned at mu_s>=%s '
                          '(spec: border nodes must be rupture-eligible)'
                          % (case, n_pinned, n_border, UNBREAKABLE))

    if fails:
        print('\nFAIL test_tpv12_tpv13_border_rupture_eligible (%d)' % len(fails))
        for f in fails:
            print(' -', f)
        return 1
    print('\nSUCCESS test_tpv12_tpv13_border_rupture_eligible: no fault-border '
          'node is pinned unbreakable in test.tpv12 or test.tpv13')
    return 0


if __name__ == '__main__':
    sys.exit(main())

#! /usr/bin/env python3
"""Per-step parity diagnostic: the Python port against Fortran's OWN internal
state, loaded from src/pydump.f90's dumps.

    run_parity.py <case_dir> [--backend numpy|jax] [--nsteps N] [--golden PATH]

This is the only thing in the tree that reads pydump_*.txt. It answers a
different question from the e2e sweep: the sweep asks "does the final answer
match?", this asks "at which step does it start to diverge?" by starting the
port from Fortran's exact pre-time-loop state instead of from a mesh this
port built itself. That is what root-caused the missing rsfNucleation branch.

ONE runner. There used to be four -- run_parity.py, run_parity_jax.py,
run_parity_rsf.py, run_parity_tp.py -- differing only in which of the six
solver modules they imported and in their argv order. friclaw and backend are
arguments now, so there is nothing left for the other three to vary.
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eqdyna import eqdyna3d                       # noqa: E402
from eqdyna import globalvar as gv                # noqa: E402
from eqdyna.pydump import load                    # noqa: E402

# frt column layout, and the fric slot each one comes from. Named rather than
# spelled as bare integers: a column read from the wrong slot is a silent
# wrong answer, not a crash.
COLUMNS = [
    ('slipS', gv.SLIP_STRIKE), ('slipD', gv.SLIP_DIP), ('slipN', gv.SLIP_NORM),
    ('srS', gv.SLIPRATE_STRIKE), ('srD', gv.SLIPRATE_DIP),
    ('peakSr', gv.SLIPRATE_MAX), ('finalSr', gv.PEAK_SLIPRATE),
    ('Tn', gv.TRACT_NORM), ('Ts', gv.TRACT_STRIKE), ('Td', gv.TRACT_DIP),
    ('vxm', gv.VEL_MASTER_X), ('vym', gv.VEL_MASTER_X + 1),
    ('vzm', gv.VEL_MASTER_X + 2),
    ('vxs', gv.VEL_SLAVE_X), ('vys', gv.VEL_SLAVE_X + 1),
    ('vzs', gv.VEL_SLAVE_X + 2),
    ('rsfstate', gv.STATE), ('thetaPc', gv.THETA_PC),
]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('case_dir')
    ap.add_argument('--backend', choices=('numpy', 'jax'), default='numpy')
    ap.add_argument('--nsteps', type=int, default=None)
    ap.add_argument('--golden', default=None,
                    help='Fortran frt text output to diff against, column by column')
    a = ap.parse_args(argv)

    t0 = time.time()
    S = load(a.case_dir)
    t1 = time.time()
    print('load %.3f s; N=%d E=%d nftnd=%d friclaw=%d'
          % (t1 - t0, S['N'], S['E'], S['nftnd'], S['friclaw']))

    out = eqdyna3d.run(S, nsteps=a.nsteps, verbose=True, backend=a.backend)
    print('solve %.3f s for %d steps (%s)'
          % (time.time() - t1, a.nsteps or S['nstep'], a.backend))

    fric = out['fric']
    coords = S['meshCoor'][S['nsmp1']]
    cols = np.column_stack(
        [coords[:, 0], coords[:, 1], coords[:, 2], out['fnft']]
        + [fric[:, slot] for _, slot in COLUMNS])
    py_out = os.path.join(a.case_dir, 'frt_py.txt')
    np.savetxt(py_out, cols, fmt='%18.7e')
    print('wrote', py_out)

    if not a.golden:
        return 0
    golden = np.loadtxt(a.golden)
    if golden.shape != cols.shape:
        print('SHAPE MISMATCH golden=%r run=%r' % (golden.shape, cols.shape))
        return 1
    labels = ['x', 'y', 'z', 'fnft'] + [n for n, _ in COLUMNS]
    print('%-22s %16s %16s' % ('col', 'max_abs_diff', 'max_rel_diff'))
    worst = 0.0
    for i, lab in enumerate(labels):
        d = np.abs(cols[:, i] - golden[:, i])
        rel = d / np.maximum(np.abs(golden[:, i]), 1e-300)
        worst = max(worst, d.max())
        print('%-22s %16.6e %16.6e' % (lab, d.max(), rel.max()))
    if np.array_equal(cols, golden):
        print('BIT-IDENTICAL to golden text output.')
    print('overall max abs diff:', worst)
    return 0


if __name__ == '__main__':
    sys.exit(main())

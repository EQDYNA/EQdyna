#! /usr/bin/env python3
"""Regression guard: the gated backend axis is exactly (fortran, python-jax)
plus the per-case python-jax-mpi opt-in -- python-numpy is in NO gate.

Owner decision 2026-09-23 ("I actually don't care numpy", "I will use Jax
anyway"): python-numpy left the everyday sweep, the release sweep and the CI
smoke. The numpy CODE stays (`python3 -m eqdyna <case> --backend numpy`), runnable
by hand; it is simply not a matrix axis, so there is no third cell state.
Mutation: add 'python-numpy' back to matrix.BACKENDS -> RED.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)


def main():
    fails = []
    try:
        from testsys import matrix
    except Exception as exc:                        # noqa: BLE001
        print('FAIL test_backend_axis: testsys.matrix failed to import (%s: %s)'
              % (type(exc).__name__, exc))
        return 1
    if tuple(matrix.BACKENDS) != ('fortran', 'python-jax', 'python-jax-mpi'):
        fails.append('matrix.BACKENDS is %r, expected exactly '
                     "('fortran', 'python-jax', 'python-jax-mpi')" % (matrix.BACKENDS,))
    runnable, unsupported = matrix.cells()
    numpy_cells = [c for c in list(runnable) + [u[:2] for u in unsupported]
                   if c[1] == 'python-numpy']
    if numpy_cells:
        fails.append('python-numpy cells present in the table: %r' % numpy_cells)
    ci_backends = sorted({b for _, b in matrix.CI_CELLS})
    if ci_backends != ['fortran', 'python-jax']:
        fails.append('matrix.CI_CELLS backends are %r, expected [fortran, python-jax]'
                     % ci_backends)
    if fails:
        print('FAIL test_backend_axis:\n  - ' + '\n  - '.join(fails))
        return 1
    print('SUCCESS test_backend_axis (%d runnable cells, backends %s, CI %s)'
          % (len(runnable), '/'.join(matrix.BACKENDS), '/'.join(ci_backends)))
    return 0


if __name__ == '__main__':
    sys.exit(main())

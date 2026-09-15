#! /usr/bin/env python3
"""
Decomposition-independent canonical form for frt output. One implementation.

WHY THIS EXISTS

`frt.txt<rank>` is written per MPI rank, so the SHAPE of a result encodes the
decomposition that produced it: test.tpv8's reference is 2 files, test.tpv29's
is 4, and the Python standalone writes exactly 1 (it is serial). Comparing them
therefore required reconciliation at every call site, and `check.test.py`
hardcoded `['frt.txt0', 'frt.txt1', 'frt.txt2', 'frt.txt3']` -- a gate whose
file list is a statement about npx*npy*npz.

The canonical form makes a result a statement about the PHYSICS instead: one
row per fault node, ordered by position. Then a 4-rank Fortran run, a serial
Fortran run, and a serial Python run are all compared the same way, against the
same artifact, with the same tolerance -- which is what makes backends
comparable quantitatively rather than each tier having its own notion of
"matches".

THE TWO OPERATIONS, and why each is needed

  1. DEDUPE by rounded (x, y, z). A node on a partition boundary is written by
     EVERY rank that owns it, so concatenating the per-rank files double-counts
     it. Measured on the committed references: 1.6% duplicated rows for
     test.tpv8, 1.0% for test.drv.a6, 3.6% for test.tpv29. A no-op for
     genuinely single-rank output.

  2. LEXSORT by (x, y, z). Rank order is arbitrary, so row order is arbitrary.
     Without sorting, a positional diff compares unrelated fault nodes and
     reports garbage that looks like a physics failure.

COORD_DECIMALS = 6 for the dedupe key: fault node coordinates are metres, and
6 decimals is a micron. Two distinct nodes are never a micron apart (the
finest gated spacing is 50 m); the SAME node written by two ranks agrees far
inside a micron because both ranks computed it from the same grid arithmetic.
So the key cannot merge distinct nodes and cannot fail to merge duplicates.

This module does NOT compare. It only produces the canonical array, so the
tolerance policy lives in one place with the caller and not in here.

PROVENANCE OF THE COMMITTED REFERENCES

`test.reference.results/<case>/frt.canonical.txt` is the committed reference
the gate compares against. Each one was produced by this module, verbatim,
from the per-rank `frt.txt*` files that the same directory used to hold:

    python3 -m testsys.frt_canonical test.reference.results/<case>

Those per-rank files were deleted in the same change, deliberately: keeping
both would have stored the same numbers twice and left the gate a choice of
two references. The transform is lossless here and that was checked rather
than assumed -- across all eight cases, 358 duplicate groups, the duplicate
rows agreed to 0.000e+00 in every one of the 22 columns (canonicalize() now
refuses to proceed if they ever do not). The row counts before and after were
tpv8 1922->1891, tpv10 1922->1891, tpv104 2738->2701, tpv1053d 4050->4005,
drv.a6 5202->5151, tpv29 3444->3321, meng2023a 672->651, meng2023cb 672->651.
Regenerating a reference is still a deliberate, reviewed commit (rule 7); the
git history holds the per-rank originals.
"""
import glob
import os

import numpy as np

COORD_DECIMALS = 6
FRT_COLUMNS = 22


def frt_rank_files(case_dir, pattern='frt.txt*'):
    """Every per-rank frt file in `case_dir`, sorted by rank number.

    Sorted NUMERICALLY, not lexically: with 10+ ranks a lexical sort gives
    frt.txt0, frt.txt1, frt.txt10, frt.txt2 ... Row order does not survive
    canonicalisation anyway, but a deterministic input order keeps the
    concatenation byte-reproducible, which is what makes the committed
    canonical file verifiable (rule 4).
    """
    paths = glob.glob(os.path.join(case_dir, pattern))

    def rank_of(p):
        tail = os.path.basename(p).split('frt.txt')[-1]
        return int(tail) if tail.isdigit() else -1

    return sorted(paths, key=rank_of)


def load_frt(paths):
    """Concatenate the given frt files into one (N, 22) array.

    Raises on an empty list rather than returning an empty array: "no files
    found" and "a run that produced no fault nodes" are different failures and
    must not look alike (rule 2).
    """
    if not paths:
        raise ValueError('load_frt: no frt files given -- refusing to return '
                         'an empty array, which would compare equal to another '
                         'empty array and report a false pass')
    blocks = []
    for p in paths:
        a = np.loadtxt(p)
        if a.ndim == 1:
            a = a.reshape(1, -1)
        if a.shape[1] != FRT_COLUMNS:
            raise ValueError('load_frt: %s has %d columns, expected %d'
                             % (p, a.shape[1], FRT_COLUMNS))
        blocks.append(a)
    return np.vstack(blocks)


def canonicalize(arr):
    """Dedupe by rounded (x,y,z), then lexsort by (x,y,z). Returns a new array.

    Deterministic: np.unique on the rounded key returns the FIRST occurrence
    index, and the subsequent lexsort fixes the order completely, so the output
    does not depend on how many ranks wrote the input or in what order they
    were concatenated.

    Dropping a duplicate is only lossless if the duplicates AGREE, so that is
    checked, not assumed: if two ranks wrote different values for the same
    fault node, "keep the first" would be a silent choice between two answers
    (rule 2), and it raises instead. Measured on all eight committed
    references (358 duplicate groups, every one of the 22 columns): max
    intra-duplicate spread 0.000e+00 -- fault-node quantities are synchronised
    across ranks, so exact equality is the right test and this cannot fire on
    ordinary floating-point noise.
    """
    if arr.ndim != 2 or arr.shape[1] != FRT_COLUMNS:
        raise ValueError('canonicalize: expected (N, %d), got %r'
                         % (FRT_COLUMNS, arr.shape))
    key = np.round(arr[:, :3], COORD_DECIMALS)
    _, first, inverse, counts = np.unique(key, axis=0, return_index=True,
                                          return_inverse=True,
                                          return_counts=True)
    inverse = inverse.reshape(-1)
    for group in np.flatnonzero(counts > 1):
        rows = arr[inverse == group]
        spread = np.abs(rows.max(axis=0) - rows.min(axis=0))
        if spread.max() > 0.0:
            col = int(spread.argmax())
            raise ValueError(
                'canonicalize: the fault node at (%.6f, %.6f, %.6f) was written '
                '%d times with DIFFERENT values -- max spread %.6e at column '
                '%d. Deduping would silently keep one of them. Rank-shared '
                'fault-node quantities are supposed to be identical; this is a '
                'finding, not a tolerance question.'
                % (rows[0, 0], rows[0, 1], rows[0, 2], rows.shape[0],
                   spread.max(), col))
    kept = arr[np.sort(first)]
    order = np.lexsort((kept[:, 2], kept[:, 1], kept[:, 0]))
    return kept[order]


def canonical_from_case(case_dir, pattern='frt.txt*'):
    """The canonical array for a case directory, whatever its rank count."""
    return canonicalize(load_frt(frt_rank_files(case_dir, pattern)))


def write_canonical(arr, path):
    """Write the canonical array in the same E18.7E4 format frt output uses.

    Same format as the per-rank files on purpose: the canonical file is then
    readable by every existing tool that already parses frt output, and a
    human diffing it against a rank file sees aligned columns.
    """
    with open(path, 'w') as fh:
        for row in arr:
            fh.write(''.join('%19.7E' % v for v in row).replace('E+0', 'E+00')
                     .replace('E-0', 'E-00') + '\n')


def align(a, b, coord_tol=1e-9):
    """Canonicalise both, then assert they describe the same node set.

    Returns (a_canon, b_canon) ready for an elementwise comparison.

    coord_tol is deliberately tiny and CANNOT mask a misalignment. For a planar
    fault both sides' coordinates are exact grid multiples and agree at 0.0.
    For a dipping or rough fault they are the OUTPUT of insertFaultInterface's
    y-blend, computed independently by Fortran and by the Python port, so a few
    ULPs (~1e-12) is expected. A genuine node-to-node misalignment shows a diff
    on the order of the grid spacing -- metres, nine orders of magnitude above
    this floor.
    """
    ac, bc = canonicalize(a), canonicalize(b)
    if ac.shape != bc.shape:
        raise ValueError(
            'align: node counts differ after canonicalisation -- %r vs %r. '
            'That is a different fault discretisation, not a tolerance '
            'question.' % (ac.shape, bc.shape))
    dev = np.abs(ac[:, :3] - bc[:, :3]).max()
    if dev > coord_tol:
        raise ValueError(
            'align: canonicalised rows do not describe the same nodes -- max '
            'coordinate difference %.3e exceeds %.3e. A real misalignment is '
            'metres; this is either a different mesh or a bug in '
            'canonicalisation.' % (dev, coord_tol))
    return ac, bc


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(
        description='Write the decomposition-independent canonical frt file '
                    'for a case directory.')
    ap.add_argument('case_dir')
    ap.add_argument('-o', '--out', default=None,
                    help='output path (default: <case_dir>/frt.canonical.txt)')
    args = ap.parse_args()
    files = frt_rank_files(args.case_dir)
    arr = canonicalize(load_frt(files))
    out = args.out or os.path.join(args.case_dir, 'frt.canonical.txt')
    write_canonical(arr, out)
    print('%d rank file(s) -> %d canonical rows -> %s'
          % (len(files), arr.shape[0], out))

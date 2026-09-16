"""Unit cover for the element-force SCATTER invariants assembleGlobalKU.py's
memory footprint now rests on (2026-09-15).

`elastic_step` used to stage the WHOLE element-force function space -- one
concatenated index array and one value buffer, 148 entries per element, 874 MB
+ 874 MB on test.tpv104 -- so that a single `np.bincount` could accumulate it.
It now accumulates block by block with `np.add.at` straight into `force`, which
removes both arrays (-1.56 GB measured on test.tpv104) and is also slightly
faster (the scatter alone: 91.5 ms vs 136.8 ms per step on test.tpv8's real
index arrays).

That substitution is only legal because of an ORDERING property, and this file
pins it. Neither testsys/parity nor testsys/accept can: both gate on a
tolerance against a reference, and the failure mode here is a last-bit
reduction-order change -- the port's output would move by ~1e-13 relative and
stay green everywhere while no longer being the same computation. The
end-to-end digest checks in the session that made the change catch it once; a
test catches it every time.

Rule 2 note: nothing here is a tolerance. Every assertion is on BYTES.
"""
import os
import sys

import numpy as np
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO_ROOT, 'src', 'python'))


def _blocks(seed, nblk=6, nrow=97, ntarget=41):
    """An index/value pattern with the two features that make the real scatter
    hard: heavy DUPLICATION within a block (nrow*8 entries into ntarget
    targets, so ~19 contributions per target) and REPEATED index arrays across
    blocks (the four hourglass modes scatter through the same three index
    arrays). Values span several magnitudes so that a reordered sum cannot
    coincidentally round the same way."""
    rng = np.random.default_rng(seed)
    idx = [rng.integers(0, ntarget, size=nrow * 8) for _ in range(3)]
    idx = (idx * ((nblk + 2) // 3))[:nblk]
    vals = [(rng.random(nrow * 8) - 0.5) * 10.0 ** rng.integers(-6, 7, size=nrow * 8)
            for _ in range(nblk)]
    return idx, vals, ntarget


@pytest.mark.parametrize('seed', [0, 1, 7])
def test_block_add_at_is_bitwise_bincount(seed):
    """THE invariant the footprint win rests on: `np.add.at` once per block, in
    block order, is BIT for BIT `np.bincount` over the concatenation of those
    blocks. Both accumulate every target's contributions in the same sequence
    (block order, then within-block order), and that is the whole argument --
    if a future NumPy changed either one's traversal, the port would silently
    stop reproducing its own reference output."""
    idx, vals, n = _blocks(seed)

    got = np.zeros(n)
    for i, v in zip(idx, vals):
        np.add.at(got, i, v)

    want = np.bincount(np.concatenate(idx), weights=np.concatenate(vals), minlength=n)

    assert got.tobytes() == want.tobytes(), \
        'np.add.at per block no longer accumulates in np.bincount order -- the ' \
        'in-place block scatter is NOT a bit-identical substitution on this NumPy'


@pytest.mark.parametrize('seed', [0, 1, 7])
def test_per_block_bincount_is_NOT_interchangeable(seed):
    """The near-miss that must never be "optimised" back in: accumulating each
    block with its OWN `np.bincount` and adding the results regroups the sum
    ((a+b)+(c+d) instead of ((a+b)+c)+d) and does NOT reproduce the reference
    bytes. Recorded as a test rather than a comment because it is the obvious
    way to write this, it is faster-looking, and it is wrong."""
    idx, vals, n = _blocks(seed)

    per_block = np.zeros(n)
    for i, v in zip(idx, vals):
        per_block += np.bincount(i, weights=v, minlength=n)

    want = np.bincount(np.concatenate(idx), weights=np.concatenate(vals), minlength=n)

    assert not np.array_equal(per_block, want), \
        'per-block bincount happened to match here; the ordering argument in ' \
        'assembleGlobalKU.build() needs re-deriving if this is no longer a real difference'


@pytest.mark.parametrize('seed', [0, 3])
def test_incremental_group_sum_matches_out_form(seed):
    """The PML group sums (f1+f2+f3v+f10 and friends) are now accumulated into
    a persistent buffer as the twelve f-blocks are built and scattered, instead
    of being computed from twelve simultaneously-live blocks. `g = f1; g += f2;
    g += f3; g += f4` is the same left-associative sum as `np.add(f1,f2,out=g);
    g += f3; g += f4`, so this is exact -- pinned because it is what lets the
    PML staging be 4 blocks instead of 15."""
    rng = np.random.default_rng(seed)
    f = [(rng.random((53, 8)) - 0.5) * 10.0 ** rng.integers(-4, 5, size=(53, 8))
         for _ in range(4)]

    incr = f[0].copy()
    for k in (1, 2, 3):
        incr += f[k]

    out = np.empty_like(f[0])
    np.add(f[0], f[1], out=out); out += f[2]; out += f[3]

    assert incr.tobytes() == out.tobytes()


def test_build_exposes_no_staged_function_space():
    """A guard against the staging coming back by accident: `inv` must not
    carry a concatenated scatter index or a whole-function-space value buffer.
    Those two keys were 1.75 GB on test.tpv104; a future refactor that
    reintroduces them would be a silent 2x on peak RSS for large meshes, and
    nothing else in the suite measures memory."""
    from eqdyna import assembleGlobalKU
    src = open(assembleGlobalKU.__file__).read()
    for gone in ("scat_idx=", "scat_val=", "scat_off="):
        assert gone not in src, \
            'assembleGlobalKU.build() is staging the element-force function space again ' \
            '(%s) -- see its in-place block-scatter note for why that costs 1.75 GB' % gone

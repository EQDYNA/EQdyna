"""Guards for the float64 identities the NumPy solver's speed depends on.

python/eqdyna/assembleGlobalKU.py's `elastic_step` was rewritten (2026-09-14) to
cut its per-step cost roughly in half. Every rewrite was chosen so the output
stays BIT-IDENTICAL, not merely close: the standalone port is gated against
Fortran at PER_CASE_ABS_BOUND['test.tpv8'] == 1e-8 while the observed diff is
8.23e-11, and a perf change has no business spending that margin.

Two of those rewrites rest on properties of NumPy itself rather than on
algebra, so they are asserted here instead of being left for the end-to-end
parity tier to notice:

  1. The hourglass contraction `np.einsum('ei,eij->ej', ...)` reduces its eight
     terms in the SAME order as the `(phi[:,m,:,None] * dl).sum(axis=1)` form it
     replaced. This is a statement about NumPy's inner loops, and it is NOT a
     given: the obvious-looking alternative, a stacked `phi @ dl` matmul, is
     measurably faster still and looked identical on a spot-check of 2.8 M real
     values -- but running it end-to-end diverged at step 5 (8.13e-20 on
     near-cancelling force entries, numpy 2.2.6, almost certainly FMA
     contraction in matmul's inner loop). That is exactly the failure mode a
     spot-check misses and a released NumPy upgrade could introduce here, so the
     surviving claim gets a test.

  2. Folding a sign into a hoisted array -- `(-phi) * r` instead of
     `-(phi * r)` -- is exact in IEEE 754, which is what lets the hourglass
     force be written straight into the scatter buffer with `out=`.

A failure here does not mean the solver is wrong; it means a NumPy upgrade has
changed a reduction order, and assembleGlobalKU.py's bit-identity claims (and the
comments asserting them) must be re-verified before shipping.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    'python'))
from eqdyna import assembleGlobalKU  # noqa: E402


def _operands(seed, n=20000):
    """Random operands PLUS a deliberately near-cancelling block.

    Uniform random data hides reduction-order differences: they only surface
    when the running sum loses most of its significant digits. The second half
    of `d` is built so the eight products very nearly cancel, which is the
    regime where the real solver diverged.
    """
    rng = np.random.default_rng(seed)
    p = rng.standard_normal((n, 8))
    d = rng.standard_normal((n, 8, 3))
    half = n // 2
    d[half:] = (rng.standard_normal((n - half, 8, 3)) * 1e-14
                + np.array([1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0])[None, :, None])
    return p, d


@pytest.mark.parametrize('seed', [0, 1, 2, 3, 4])
def test_hourglass_einsum_matches_broadcast_reduce_bitwise(seed):
    """einsum('ei,eij->ej') == (p[:,:,None]*d).sum(axis=1), bit for bit."""
    p, d = _operands(seed)
    reference = (p[:, :, None] * d).sum(axis=1)
    fast = np.einsum('ei,eij->ej', p, d)
    ndiff = int((reference != fast).sum())
    assert ndiff == 0, (
        'np.einsum no longer reduces in the same order as the broadcast-multiply'
        '/axis-1 sum it replaced in assembleGlobalKU.assembleGlobalKU: %d of %d values '
        'differ, max abs %.3e. The hourglass contraction is no longer '
        'bit-identical -- re-verify the port against Fortran before shipping.'
        % (ndiff, reference.size, np.abs(reference - fast).max()))


def test_negating_a_factor_is_exact():
    """(-a)*b == -(a*b) for every finite double: only the sign bit moves."""
    rng = np.random.default_rng(7)
    a = rng.standard_normal(200000) * 10.0 ** rng.integers(-100, 100, 200000)
    b = rng.standard_normal(200000) * 10.0 ** rng.integers(-100, 100, 200000)
    assert np.array_equal((-a) * b, -(a * b))


@pytest.mark.parametrize('seed', [0, 1, 2])
def test_c_contraction_is_layout_independent(seed):
    """`_c` is unchanged by making its operands contiguous.

    build() now hands `_c` contiguous copies of dNx/dNy/dNz and per-component
    velocity gathers instead of strided views of an (E,8,3) block. That is a
    layout change only: `dN * v` allocates a contiguous temporary either way, so
    `.sum(axis=1)` sees identical bytes and reduces them identically.
    """
    rng = np.random.default_rng(seed)
    blk_dn = rng.standard_normal((30000, 8, 3))
    blk_v = rng.standard_normal((30000, 8, 3))
    for k in range(3):
        strided = assembleGlobalKU._c(np, blk_dn[:, :, k], blk_v[:, :, k])
        contiguous = assembleGlobalKU._c(np, np.ascontiguousarray(blk_dn[:, :, k]),
                                      np.ascontiguousarray(blk_v[:, :, k]))
        assert np.array_equal(strided, contiguous)


def test_c_is_not_interchangeable_with_einsum():
    """`_c`'s docstring says einsum is NOT a drop-in for it. Pin that down.

    `_c` reduces a contiguous (E,8) temporary, so NumPy uses pairwise summation;
    einsum accumulates sequentially. They disagree in the last bits. This test
    exists so that anyone "optimising" `_c` into the einsum spelling -- the one
    substitution that looks obviously right and is 2.7x faster in isolation --
    finds out here rather than from a parity run.

    It asserts that a difference EXISTS on a cancellation-heavy input. If NumPy
    ever makes the two agree, this test fails, and the right response is to
    re-read `_c`'s docstring and re-measure -- not to delete the test.
    """
    rng = np.random.default_rng(11)
    n = 200000
    dN = rng.standard_normal((n, 8))
    v = rng.standard_normal((n, 8)) * 1e-13 + np.tile([1.0, -1.0], 4)[None, :]
    assert not np.array_equal(assembleGlobalKU._c(np, dN, v), np.einsum('ei,ei->e', dN, v))

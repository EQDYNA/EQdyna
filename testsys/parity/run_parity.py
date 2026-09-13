#! /usr/bin/env python3
"""
`parity` testsys tier (pathway_forward.md item 14). Runs python/eqdyna's
NumPy and (if importable) JAX tpv8 ports against the fixtures
testsys/parity/make_fixtures.py generates, and FAILS (non-zero exit) if
any output column exceeds the roundoff/ASCII-quantization envelope
established in python/README-parity.md's Update 1-3 (tpv8, friclaw=1,
114-step full run vs the golden `frt.txt`).

Does NOT silently pass if the Fortran binary/fixtures are missing --
raises loudly, per PROJECT_RULES rule 2/3.

Thresholds (see THRESHOLDS below) are 2-5x the actual observed max abs
diff from three independent runs logged in README-parity.md Updates 1-3
and 6 (49.54 Pa max on Tn/Ts, ~5e-7 s max on fnft) -- generous enough to
absorb ordinary machine-to-machine floating-point reduction-order noise,
tight enough that a real algorithmic regression (which historically has
shown up as >1 order of magnitude, not a 2-5x wobble) still fails the
gate. Columns the Fortran never writes for friclaw=1 (fric 31-36/47/20/23)
are gated at exactly 0 -- any non-zero value there means something now
writes to a column that structurally should not be touched for this
friction law, which is exactly the kind of silent regression this tier
exists to catch.
"""
import os
import sys

import numpy as np

TESTSYS = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(TESTSYS))
PYTHON_PKG = os.path.join(REPO_ROOT, 'python')
sys.path.insert(0, PYTHON_PKG)
sys.path.insert(0, os.path.join(PYTHON_PKG, 'eqdyna'))

FIXTURE_CASE = os.path.join(TESTSYS, 'fixtures', 'test_tpv8_serial')

# column label -> max abs diff allowed (see module docstring for provenance)
THRESHOLDS = {
    'x': 0.0, 'y': 0.0, 'z': 0.0,                 # coordinates: must be bit-identical
    'fnft(rupt.time)': 2e-6,
    'slipS': 2e-6, 'slipD': 5e-8, 'slipN': 1e-14,
    'srS': 2e-6, 'srD': 5e-8, 'peakSr': 2e-6,
    'finalSr(47,unused)': 0.0,                    # never written for friclaw=1
    'Tn(78)': 200.0, 'Ts(79)': 200.0, 'Td(80)': 2.0,
    'vxm(31)': 0.0, 'vym(32)': 0.0, 'vzm(33)': 0.0,   # never written for friclaw=1
    'vxs(34)': 0.0, 'vys(35)': 0.0, 'vzs(36)': 0.0,
    'rsfstate(20,unused)': 0.0, 'thetaPc(23,unused)': 0.0,
}
LABELS = list(THRESHOLDS.keys())


def _frt_columns(out, S):
    fric = out['fric']
    nsmp1 = S['nsmp1']
    coords = S['meshCoor'][nsmp1]
    return np.column_stack([
        coords[:, 0], coords[:, 1], coords[:, 2], out['fnft'],
        fric[:, 70], fric[:, 71], fric[:, 72], fric[:, 73], fric[:, 74], fric[:, 75],
        fric[:, 46], fric[:, 77], fric[:, 78], fric[:, 79],
        fric[:, 30], fric[:, 31], fric[:, 32], fric[:, 33], fric[:, 34], fric[:, 35],
        fric[:, 19], fric[:, 22],
    ])


def check_engine(name, run_fn, S, golden, nsteps):
    print(f'-- parity: {name} --')
    with np.errstate(divide='ignore', invalid='ignore'):
        out = run_fn(S, nsteps=nsteps, verbose=False)
    cols = _frt_columns(out, S)
    if cols.shape != golden.shape:
        print(f'FAIL {name}: shape mismatch {cols.shape} vs golden {golden.shape}')
        return False
    ok = True
    for i, label in enumerate(LABELS):
        d = float(np.abs(cols[:, i] - golden[:, i]).max())
        limit = THRESHOLDS[label]
        passed = d <= limit
        ok = ok and passed
        flag = 'ok' if passed else 'EXCEEDS THRESHOLD'
        print(f'  {label:22s} max_abs_diff={d:.6e}  limit={limit:.6e}  {flag}')
    print(f'{"SUCCESS" if ok else "FAIL"} parity:{name}')
    return ok


def main():
    if not os.path.isdir(FIXTURE_CASE):
        raise SystemExit(
            f'FAIL: fixtures missing at {FIXTURE_CASE}. '
            'Run: python3 testsys/parity/make_fixtures.py')
    golden_path = os.path.join(FIXTURE_CASE, 'frt.txt0')
    if not os.path.exists(golden_path):
        raise SystemExit(f'FAIL: golden oracle missing: {golden_path}')

    from eqdyna.port import load  # noqa: E402 -- import after sys.path setup
    S = load(FIXTURE_CASE)
    golden = np.loadtxt(golden_path)
    nsteps = S['nstep']

    results = {}
    from eqdyna import port
    results['numpy'] = check_engine('numpy', port.run, S, golden, nsteps)

    try:
        from eqdyna import port_jax
    except ImportError as e:
        print(f'SKIP parity:jax -- jax not importable ({e}); numpy result alone gates this tier.')
    else:
        results['jax'] = check_engine('jax', port_jax.run, S, golden, nsteps)

    overall = 0 if all(results.values()) else 1
    print(f'\n{"SUCCESS" if overall == 0 else "FAIL"} parity (engines checked: {list(results)})')
    return overall


if __name__ == '__main__':
    sys.exit(main())

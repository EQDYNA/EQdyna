#! /usr/bin/env python3
"""Row 166 probe (rule 4a), by hand: does the row-165 locked-node variant ever
put the slip-rate magnitude that faulting divides by at 0 or at a subnormal?

Suspect (a): solveRSF's unit direction srS/|sr| (faulting.f90:296,303,
faulting.py trialS/Ts) is 0/0 at V=0.  Suspect (b): subnormal flush (FTZ)
changes the amplified roundoff.  Both need |sr| to reach 0 or < 2.3e-308.
This runs the variant on the numpy backend (serial, no FTZ, no jit hiding the
arrays) and records, on every call of fric.rate_state_ageing_law, the minimum of
the pre-loop v_trial = sqrt((srS+vini_s)^2+(srD+vini_d)^2) over ALL fault nodes.

MEASURED 2026-11-06, 240 steps (5 s): 5040 calls, exact zeros 0, subnormals 0,
count below 1e-150 0, minimum 1.07e-12 (the background creep vini keeps even
the slip_rate=1e-30 locked rows at >= 1e-12).  So neither the 0/0 branch nor
a subnormal is ever reached: both suspects FALSIFIED on this observable.  jax
is always float64 (backend.array_module sets jax_enable_x64), so "x64 not run"
was moot.

Usage: test_row166_locked_v_floor.py [nsteps=2] [--floor F]   (default F=1e-150;
`--floor 1e-9` must FAIL: the floor is 1.07e-12, so the check can go red.)
Excluded from the sweep (by hand; runs the numpy solver, ~2 s/step).
"""
import os
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import row165_locked_variant as R  # noqa: E402

DT = 5.0 / 240


def main(argv):
    floor = 1e-150
    if '--floor' in argv:
        i = argv.index('--floor')
        floor = float(argv[i + 1])
        del argv[i:i + 2]
    nt = int(argv[0]) if argv else 2
    sys.path.insert(0, os.path.join(R.REPO, 'src', 'python'))
    mins = []
    from eqdyna import fric as F
    orig = F.rate_state_ageing_law

    def wrapped(xp, v, *a, **k):
        mins.append(float(np.asarray(v).min()))
        return orig(xp, v, *a, **k)
    F.rate_state_ageing_law = wrapped
    with tempfile.TemporaryDirectory() as root:
        fx = os.path.join(root, 'variant.nc')
        R.make_variant(fx)
        cd = R.stage(root, 'p', fx, (nt - 0.5) * DT)
        with open(os.path.join(cd, 'user_defined_params.py'), 'a') as f:
            f.write('\npar.nx = 1\npar.ny = 1\npar.nz = 1\n')
        env = dict(R._env(), PYTHONPATH=os.path.join(R.REPO, 'src', 'python'))
        subprocess.run([sys.executable, 'case.setup'], cwd=cd, env=env, check=True,
                       stdout=subprocess.DEVNULL)
        os.environ['JAX_PLATFORMS'] = 'cpu'
        sys.argv = ['eqdyna', cd, '--backend', 'numpy', '--device', 'cpu']
        from eqdyna.eqdyna3d import main as run_main
        run_main()
    if not mins:
        print('FAIL: rate_state_ageing_law was never called (probe measured nothing)')
        return 1
    m = min(mins)
    print('calls %d, min |sr| over all fault nodes %.3e (floor %.1e)' % (len(mins), m, floor))
    if not m > floor:
        print('FAIL: |sr| reached %.3e <= %.1e: suspect (a)/(b) is live' % (m, floor))
        return 1
    print('SUCCESS: |sr| never reaches 0 or a subnormal; row-166 suspects stay falsified')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

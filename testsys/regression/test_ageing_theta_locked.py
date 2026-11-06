#! /usr/bin/env python3
"""
Regression guard (PROJECT_RULES.md rule 10; board row 163): the ageing-law
state update must stay finite and positive at locked-node slip rates, in BOTH
implementations.

The old update  theta = L/V + (theta - L/V)*exp(-V*dt/L)  cancels
catastrophically: at V ~ 1e-30 (the 1e-30 the code adds to V only in the xmu
line), L/V ~ 1e29, exp(-V dt/L) = 1 exactly and theta collapses to exactly 0;
at V == 0, L/0 = inf and inf - inf = NaN, which the next call's
dlog(V0*theta/L) and the Newton update spread to the slip rate.

Checks (behaviour only, rule 10a):
  1. python NewtonRaphson (faulting.NewtonRaphson), a locked node at v_trial=0:
     slip rate and state stay finite for 3 steps;
  2. python NewtonRaphson, a node at 1e-30 with no driving traction (the slip
     rate halves each Newton iteration, 20 per step) for 60 steps, past the
     point (~54 steps) where it underflows to exactly 0: all finite and the
     state stays within 1% of theta0 + steps*dt (ageing, V ~ 0);
  3. the Fortran rate_state_ageing_law (fric.f90, compiled with globalvar.f90
     into a driver by this test) chained 3 steps at V=0 and 60 steps at 1e-30:
     theta finite, xmu finite, theta ~ theta0 + steps*dt;
  4. Fortran and python agree (rel 1e-13) with each other and with an
     independent oracle, theta*exp(-x) + dt*(-expm1(-x))/x, over V from 0 to
     1 m/s, including x straddling the 1e-3 series threshold. Both numbers
     are printed. The series threshold |x| = V*dt/L = 1e-3 is V = 0.14 m/s
     here (L=0.14, dt=1e-3); the V list straddles it on both signs, so a
     negative V (|x| large) must take the exact branch (the abs in the guard).
  5. checks 1-2 and the python half of 3-4 run on numpy AND on jax (CPU, via
     backend.array_module('jax'), the gated backend); jax == numpy to RTOL.
Fails (never skips) when gfortran is missing.
"""
import math
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, 'src', 'python'))
from eqdyna import backend as BK, faulting as F, fric as FR, globalvar as gv  # noqa: E402

A, B, L, F0, V0, DT = 0.01, 0.014, 0.14, 0.6, 1e-6, 1e-3
RTOL = 1e-13

DRIVER = """program drv
  use globalvar
  implicit none
  real (kind = dp) :: fricsgl(100), v, th, xmu, dxmudv, th0
  integer :: n, i, ios
  fricsgl = 0.0d0
  fricsgl(FRIC_SLOT_RSF_A) = %.17e_dp
  fricsgl(FRIC_SLOT_RSF_B) = %.17e_dp
  fricsgl(FRIC_SLOT_RSF_DC) = %.17e_dp
  fricsgl(FRIC_SLOT_RSF_R0) = %.17e_dp
  fricsgl(FRIC_SLOT_RSF_V0) = %.17e_dp
  dt = %.17e_dp
  do
    read(*,*,iostat=ios) v, th0, n
    if (ios /= 0) exit
    th = th0
    do i = 1, n
      call rate_state_ageing_law(v, th, fricsgl, xmu, dxmudv)
    enddo
    write(*,'(3ES26.17E3)') th, xmu, dxmudv
  enddo
end program drv
""" % (A, B, L, F0, V0, DT)


def build_fortran(tmp):
    if shutil.which('gfortran') is None:
        raise RuntimeError('gfortran not found')
    src = os.path.join(ROOT, 'src', 'fortran')
    drv = os.path.join(tmp, 'drv.f90')
    with open(drv, 'w') as fh:
        fh.write(DRIVER)
    exe = os.path.join(tmp, 'drv')
    subprocess.run(['gfortran', '-O3', '-ffree-line-length-none', '-o', exe,
                    os.path.join(src, 'globalvar.f90'),
                    os.path.join(src, 'fric.f90'), drv],
                   cwd=tmp, check=True, capture_output=True, timeout=300)
    return exe


def fortran(exe, cases):
    """cases: [(V, theta0, nsteps)] -> [(theta, xmu, dxmudv)]"""
    inp = ''.join('%.17e %.17e %d\n' % c for c in cases)
    p = subprocess.run([exe], input=inp, capture_output=True, text=True,
                       check=True, timeout=60)
    rows = [tuple(float(t) for t in ln.split()) for ln in p.stdout.splitlines()]
    assert len(rows) == len(cases), (len(rows), len(cases))
    return rows


def fric_arr(n=1):
    f = np.zeros((n, max(gv.RSF_A, gv.RSF_B, gv.RSF_DC, gv.RSF_R0, gv.RSF_V0,
                         gv.COHESION, gv.STATE) + 1))
    f[:, gv.RSF_A] = A; f[:, gv.RSF_B] = B; f[:, gv.RSF_DC] = L
    f[:, gv.RSF_R0] = F0; f[:, gv.RSF_V0] = V0
    return f


def py_chain(v, th0, n, xp=np):
    f = xp.asarray(fric_arr())
    th = xp.asarray(np.array([th0]))
    with np.errstate(all='ignore'):
        for _ in range(n):
            xmu, dxmudv, th = FR.rate_state_ageing_law(xp, xp.asarray(np.array([v])), th, f, DT)
    return float(th[0]), float(xmu[0]), float(dxmudv[0])


def newton_chain(v0, theta0, trial_mag, nsteps, xp=np):
    """faulting.NewtonRaphson stepped like the time loop: state0 <- state_out."""
    a = xp.asarray
    f = a(fric_arr())
    v = a(np.array([v0])); th = a(np.array([theta0])); tp = a(np.array([50e6]))
    tn = a(np.array([-50e6])); tc = a(np.array([1e-9])); tm = a(np.array([trial_mag]))
    hist = []
    with np.errstate(all='ignore'):
        for _ in range(nsteps):
            v, _, th, _ = F.NewtonRaphson(xp, 3, f, v, th, tp, tn, tm, tc, DT)
            hist.append((float(v[0]), float(th[0])))
    return hist


def oracle(v, th, dt=DT):
    x = v * dt / L
    g = 1.0 if x == 0.0 else -math.expm1(-x) / x
    return th * math.exp(-x) + dt * g


def check(exe):
    fails = []
    XPS = [('numpy', np), ('jax', BK.array_module('jax'))]

    def bad(msg):
        fails.append(msg)

    for xn, xp in XPS:
        # 1. python NewtonRaphson, locked node at v_trial == 0
        h = newton_chain(0.0, 1e5, 1e5, 3, xp)
        if not all(np.isfinite(a) and np.isfinite(b) for a, b in h):
            bad('%s NewtonRaphson, v_trial=0, locked: NaN/inf in (v, theta) %r' % (xn, h))
        # 2. python NewtonRaphson, V=1e-30, no driving traction, 60 steps
        n2, th0 = 60, 1e5
        h = newton_chain(1e-30, th0, 0.0, n2, xp)
        if not all(np.isfinite(a) and np.isfinite(b) for a, b in h):
            k = next(i for i, (a, b) in enumerate(h) if not (np.isfinite(a) and np.isfinite(b)))
            bad('%s NewtonRaphson, v_trial=1e-30, 60 steps: NaN/inf first at step %d: %r'
                % (xn, k + 1, h[k]))
        elif not all(abs(b - (th0 + (i + 1) * DT)) < 0.01 * th0 for i, (a, b) in enumerate(h)):
            bad('%s NewtonRaphson, v_trial=1e-30: state left theta0 + steps*dt by >1%%: %r'
                % (xn, h[-1]))
    # 3. Fortran and python chains on the locked inputs
    for v, th0, n in [(0.0, 1e5, 3), (1e-30, 1e5, 60), (1e-30, 1e9, 1)]:
        ft, fx, fd = fortran(exe, [(v, th0, n)])[0]
        pt, px, pd = py_chain(v, th0, n)
        jt, jx, _ = py_chain(v, th0, n, XPS[1][1])
        if not (jt == pt or abs(jt - pt) <= RTOL * abs(pt)):
            bad('jax theta %.17e != numpy %.17e at V=%g n=%d' % (jt, pt, v, n))
        want = th0 + n * DT
        print('  V=%g theta0=%g n=%d  fortran theta=%.17e xmu=%.6e  python theta=%.17e xmu=%.6e'
              % (v, th0, n, ft, fx, pt, px))
        for tag, t, x in (('fortran', ft, fx), ('python', pt, px), ('jax', jt, jx)):
            if not (math.isfinite(t) and math.isfinite(x)):
                bad('%s chain V=%g n=%d: theta/xmu not finite (%r, %r)' % (tag, v, n, t, x))
            elif abs(t - want) > 1e-6 * th0:
                bad('%s chain V=%g n=%d: theta %.9e != theta0+n*dt %.9e' % (tag, v, n, t, want))
    # 4. Fortran vs python vs oracle across V
    pos = [0.0, 1e-30, 1e-20, 1e-12, 1e-9, 1e-7, 1e-5, 1e-3, 1e-1,
           0.139, 0.1399, 0.1401, 0.141, 1.0]
    vs = pos + [-v for v in pos[1:]] + [-10.0, -1e3]
    cases = [(v, th, 1) for v in vs for th in (1e2, 1e5, 1e9)]
    frows = fortran(exe, cases)
    worst = 0.0
    for (v, th, _), (ft, _, _) in zip(cases, frows):
        pt = py_chain(v, th, 1)[0]
        jt = py_chain(v, th, 1, XPS[1][1])[0]
        ot = oracle(v, th)
        for tag, t in (('fortran', ft), ('python', pt), ('jax', jt)):
            e = abs(t - ot) / abs(ot)
            worst = max(worst, e)
            if not e <= RTOL:
                bad('%s V=%g theta=%g: %.17e vs oracle %.17e (rel %.2e)'
                    % (tag, v, th, t, ot, e))
        if abs(ft - pt) > RTOL * abs(ot):
            bad('fortran %.17e != python %.17e at V=%g theta=%g' % (ft, pt, v, th))
    print('  %d (V,theta) points, worst rel vs expm1 oracle = %.2e' % (len(cases), worst))
    return fails


def main():
    try:
        with tempfile.TemporaryDirectory() as tmp:
            fails = check(build_fortran(tmp))
    except Exception as e:  # noqa: BLE001 -- report, never skip
        fails = ['%s: %s' % (type(e).__name__, e)]
    for f in fails:
        print('  FAIL ', f)
    if not fails:
        print('  PASS  ageing-law theta finite at V=0 and 1e-30, Fortran == python == oracle')
    print('%s: %d failure(s)' % ('SUCCESS' if not fails else 'FAILURE', len(fails)))
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())

#!/usr/bin/env python3
"""ana.py <run_dir> [t=4.0] : dp000 n-stress at t, and diagTrace decomposition if present."""
import sys, os, glob, time
import numpy as np

d = sys.argv[1]
tq = float(sys.argv[2]) if len(sys.argv) > 2 else 4.0


def station(fn):
    rows = [l.split() for l in open(fn) if l.strip() and l.lstrip()[0] not in '#t']
    return np.array(rows, float)


for st in ('faultst000dp000.txt', 'faultst000dp015.txt', 'faultst000dp045.txt'):
    fn = os.path.join(d, st)
    if not os.path.exists(fn):
        continue
    a = station(fn)
    i = np.argmin(abs(a[:, 0] - tq))
    print(f"{fn} mtime={time.strftime('%Y-%m-%d %H:%M', time.localtime(os.path.getmtime(fn)))} "
          f"t={a[i,0]:.4f} vslip={a[i,4]:.4f} vrate={a[i,5]:.4f} vshear={a[i,6]:.4f} nstress={a[i,7]:.4f} "
          f"(t0 nstress={a[0,7]:.4f}, nrows={len(a)})")

fs = sorted(glob.glob(os.path.join(d, 'diagTrace.*.txt')))
rows = []
for f in fs:
    if os.path.getsize(f):
        rows.append(np.loadtxt(f, ndmin=2))
if not rows:
    sys.exit()
A = np.vstack(rows)
A[:, 2] = np.round(A[:, 2], 3) + 0.0
A[:, 1] = np.round(A[:, 1], 9)
nt, t, z, arn, mS, mM = A[:, 0], A[:, 1], A[:, 2], A[:, 3], A[:, 4], A[:, 5]
un, ud = A[:, 6:9], A[:, 9:12]
kuS, kuM, hgS, hgM = A[:, 12:15], A[:, 15:18], A[:, 18:21], A[:, 21:24]
dv, du = A[:, 24:27], A[:, 27:30]
dt = np.median(np.diff(np.unique(t)))
den = (mS + mM) * arn


def comb(fS, fM, e):
    return (mS * (fM * e).sum(1) - mM * (fS * e).sum(1)) / den / 1e6


def kin(e):
    return mS * mM * ((dv * e).sum(1) + (du * e).sum(1) / dt) / dt / den / 1e6


i0 = np.argmin(abs(np.unique(t) - tq))
tt = np.unique(t)[i0]
print(f"diagTrace: {len(A)} rows, {len(np.unique(z))} depths, dt={dt:.5f}, t={tt:.4f}")
zs = np.array(sorted(np.unique(z), reverse=True))
sd = -zs / np.sin(np.radians(60.0))          # down-dip distance of each row
prof = {}
for nm_, (fS, fM) in (('KU', (kuS, kuM)), ('HG', (hgS, hgM))):
    c = comb(fS, fM, un)
    prof[nm_] = np.array([c[np.where((z == zz) & (t == tt))[0][0]] for zz in zs])
prof['tot'] = prof['KU'] + prof['HG']
for L in (250., 500., 1000.):
    if sd[-1] < L:
        continue
    s = np.linspace(0, L, 2001)
    print(f"  mean Tn over top {L:6.0f} m down-dip (piecewise-linear nodal): " +
          "  ".join(f"{k}={np.trapezoid(np.interp(s, sd, v), s)/L:8.3f}" for k, v in prof.items()))
print("   z(m)    arn      mS        mM     | Tn: KU      HG     kin    total | Td: KU      HG     kin    total | open(m)  slip_d(m)")
for zz in sorted(np.unique(z), reverse=True):
    m = (z == zz) & (t == tt)
    k = np.where(m)[0][0]
    s = slice(k, k + 1)
    out = []
    for e in (un, ud):
        a_, b_, c_ = comb(kuS, kuM, e)[k], comb(hgS, hgM, e)[k], kin(e)[k]
        out += [a_, b_, c_, a_ + b_ + c_]
    print(f"{zz:8.1f} {arn[k]:8.1f} {mS[k]:9.3e} {mM[k]:9.3e} | "
          f"{out[0]:7.3f} {out[1]:7.3f} {out[2]:7.3f} {out[3]:7.3f} | "
          f"{out[4]:7.3f} {out[5]:7.3f} {out[6]:7.3f} {out[7]:7.3f} | "
          f"{(du[k]*un[k]).sum():9.2e} {(du[k]*ud[k]).sum():8.4f}")

#!/usr/bin/env python3
"""Williams-type corner eigenvalues for a fault meeting a free surface (plane strain).

Two elastic wedges of the same material, apex at the fault trace: wedge 1 from the
free surface (theta = 0) to the fault (theta = alpha), wedge 2 from the fault to the
other free surface (theta = pi). Airy function r^(lam+1) F(theta) in each wedge, so
stress ~ r^(lam-1): a root with 0 < Re(lam) < 1 is a stress singularity at the trace,
and a nodal traction sampled at distance ~dx scales as dx^(lam-1).
Interface: sigma_tt, sigma_rt and u_t continuous, plus either
  bonded:   u_r continuous (must give integer roots only: validation), or
  sliding:  sigma_rt = s*f*sigma_tt (f = 0 is the code's branch for a trace node in
            tension: shear held at the constant cohesion, no opening).
usage: corner_eigen.py [nu]
"""
import sys
import numpy as np

NU = float(sys.argv[1]) if len(sys.argv) > 1 else 0.25
KAP = 3.0 - 4.0 * NU


def fF(lam, th):
    return np.array([np.cos((lam + 1) * th), np.sin((lam + 1) * th), np.cos((lam - 1) * th), np.sin((lam - 1) * th)])


def fdF(lam, th):
    return np.array([-(lam + 1) * np.sin((lam + 1) * th), (lam + 1) * np.cos((lam + 1) * th),
                     -(lam - 1) * np.sin((lam - 1) * th), (lam - 1) * np.cos((lam - 1) * th)])


def fut(lam, th):
    return np.array([(lam + 1) * np.sin((lam + 1) * th), -(lam + 1) * np.cos((lam + 1) * th),
                     (KAP + lam) * np.sin((lam - 1) * th), -(KAP + lam) * np.cos((lam - 1) * th)])


def fur(lam, th):
    return np.array([-(lam + 1) * np.cos((lam + 1) * th), -(lam + 1) * np.sin((lam + 1) * th),
                     (KAP - lam) * np.cos((lam - 1) * th), (KAP - lam) * np.sin((lam - 1) * th)])


def det(lam, alpha, mode, f, s):
    lam = complex(lam)
    z = np.zeros(4, dtype=complex)
    m = np.zeros((8, 8), dtype=complex)
    m[0] = np.concatenate([fF(lam, 0.0), z])
    m[1] = np.concatenate([fdF(lam, 0.0), z])
    m[2] = np.concatenate([z, fF(lam, np.pi)])
    m[3] = np.concatenate([z, fdF(lam, np.pi)])
    m[4] = np.concatenate([fF(lam, alpha), -fF(lam, alpha)])
    m[5] = np.concatenate([fdF(lam, alpha), -fdF(lam, alpha)])
    m[6] = np.concatenate([fut(lam, alpha), -fut(lam, alpha)])
    if mode == 'bonded':
        m[7] = np.concatenate([fur(lam, alpha), -fur(lam, alpha)])
    else:
        m[7] = np.concatenate([fdF(lam, alpha) + s * f * (lam + 1) * fF(lam, alpha), z])
    return np.linalg.det(m)


def roots(alpha, mode, f=0.0, s=1.0):
    found = []
    for re0 in np.arange(0.05, 2.05, 0.05):
        for im0 in (0.0, 0.1, 0.3, 0.6):
            x0 = complex(re0, im0)
            x1 = x0 + 0.01
            d0 = det(x0, alpha, mode, f, s)
            d1 = det(x1, alpha, mode, f, s)
            ok = False
            for _ in range(80):
                if d1 == d0:
                    break
                x2 = x1 - d1 * (x1 - x0) / (d1 - d0)
                x0, d0 = x1, d1
                x1 = x2
                d1 = det(x1, alpha, mode, f, s)
                if abs(x1 - x0) < 1e-12:
                    ok = True
                    break
            if not ok or not (0.02 < x1.real < 2.05) or abs(x1.imag) > 1.5:
                continue
            x1 = complex(x1.real, abs(x1.imag))
            if abs(x1.imag) < 1e-8:
                x1 = complex(x1.real, 0.0)
            if all(abs(x1 - q) > 1e-5 for q in found):
                found.append(x1)
    return sorted(found, key=lambda q: (q.real, q.imag))


def show(tag, alpha_deg, mode, f=0.0, s=1.0):
    r = roots(np.radians(alpha_deg), mode, f, s)
    txt = ', '.join('%.4f%+.4fi' % (q.real, q.imag) if q.imag else '%.4f' % q.real for q in r)
    print('%-52s lam: %s' % (tag, txt))
    return r


if __name__ == '__main__':
    print('nu = %.3f, kappa = %.3f' % (NU, KAP))
    show('validation: bonded, alpha 60 (half-plane)', 60, 'bonded')
    show('validation: bonded, alpha 90 (half-plane)', 90, 'bonded')
    for a in (60, 90, 120):
        show('sliding f=0 (tension branch), alpha %d' % a, a, 'slide', 0.0, 1.0)
        for f in (0.1, 0.7):
            for s in (1.0, -1.0):
                show('sliding f=%.1f s=%+d, alpha %d' % (f, s, a), a, 'slide', f, s)


def smin(lam, alpha, f, s):
    """smallest/largest singular value of the sliding-interface matrix (0 at a root)."""
    lam = complex(lam)
    z = np.zeros(4, dtype=complex)
    m = np.array([np.concatenate([fF(lam, 0.0), z]), np.concatenate([fdF(lam, 0.0), z]),
                  np.concatenate([z, fF(lam, np.pi)]), np.concatenate([z, fdF(lam, np.pi)]),
                  np.concatenate([fF(lam, alpha), -fF(lam, alpha)]),
                  np.concatenate([fdF(lam, alpha), -fdF(lam, alpha)]),
                  np.concatenate([fut(lam, alpha), -fut(lam, alpha)]),
                  np.concatenate([fdF(lam, alpha) + s * f * (lam + 1) * fF(lam, alpha), z])])
    sv = np.linalg.svd(m, compute_uv=False)
    return sv[-1] / sv[0]


def scan(alpha_deg, f, s, re_lo, re_hi):
    """Completeness check independent of the secant search: grid minima of smin in a strip."""
    al = np.radians(alpha_deg)
    res = np.arange(re_lo, re_hi + 1e-9, 0.005)
    ims = np.arange(0.0, 1.5 + 1e-9, 0.005)
    g = np.array([[smin(complex(a, b), al, f, s) for b in ims] for a in res])
    hits = []
    for i in range(1, len(res) - 1):
        for j in range(len(ims)):
            v = g[i, j]
            nb = [g[i - 1, j], g[i + 1, j]]
            if j > 0:
                nb.append(g[i, j - 1])
            if j < len(ims) - 1:
                nb.append(g[i, j + 1])
            if v < 2e-3 and v <= min(nb):
                hits.append('%.3f%+.3fi' % (res[i], ims[j]))
    print('scan alpha %3d f=%.1f s=%+d Re in [%.2f, %.2f]: min smin %.2e; grid minima below 2e-3: %s'
          % (alpha_deg, f, s, re_lo, re_hi, g.min(), ', '.join(hits) if hits else 'none'))


if __name__ == '__main__':
    print('completeness scans (Re < 1 would be a stress singularity)')
    for f, s in ((0.0, 1.0), (0.1, 1.0), (0.1, -1.0)):
        scan(60, f, s, 0.05, 0.98)
        scan(60, f, s, 1.02, 1.70)
    scan(90, 0.1, 1.0, 0.05, 0.98)

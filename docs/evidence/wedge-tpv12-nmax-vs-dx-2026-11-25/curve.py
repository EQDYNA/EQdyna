"""n_max vs dx on the TPV12 trace node (faultst000dp000, col 8 = normal stress, + = tension).
usage: curve.py LABEL=path/to/faultst000dp000.txt:DX ...   (prints a table and the fits)"""
import sys, numpy as np
rows = []
for a in sys.argv[1:]:
    lab, rest = a.split("=", 1); path, dx = rest.rsplit(":", 1)
    d = []
    for l in open(path):
        try: d.append([float(x) for x in l.split()])
        except ValueError: pass          # header / comment lines
    d = np.array(d)
    t, n = d[:, 0], d[:, 7]; m = t <= 4.5; t, n, dd = t[m], n[m], d[m]
    i = n.argmax(); sr = np.hypot(dd[:, 2], dd[:, 5]); tau = np.hypot(dd[:, 3], dd[:, 6])
    rows.append((lab, float(dx), len(t), n[i], t[i], sr[i], sr.max(), tau[i], n[np.abs(t-4).argmin()]))
    print(f"{lab:6s} dx {float(dx):5.0f} samples {len(t):4d} n_max {n[i]:+8.4f} MPa t {t[i]:.3f}  sr@peak {sr[i]:.3f} sr_max {sr.max():.3f} m/s  |tau|@peak {tau[i]:.3f} MPa  n(4s) {rows[-1][-1]:+.3f}")
w = [r for r in rows if r[0].startswith("wedge")]
if len(w) >= 3:
    dx = np.array([r[1] for r in w]); y = np.array([r[3] for r in w]); x = np.log(500.0/dx)
    A = np.vstack([x, np.ones_like(x)]).T; c, *_ = np.linalg.lstsq(A, y, rcond=None); res = y - A@c
    print(f"log fit n_max = {c[1]:.3f} + {c[0]:.3f}*ln(500/dx): residuals {np.round(res,3)} max|res| {abs(res).max():.3f} MPa")
    for p in (0.5, 1.0):
        x2 = (500.0/dx)**p - 1; A2 = np.vstack([x2, np.ones_like(x2)]).T; c2, *_ = np.linalg.lstsq(A2, y, rcond=None); r2 = y - A2@c2
        print(f"power p={p}: n_max = {c2[1]:.3f} + {c2[0]:.3f}*((500/dx)^{p}-1): max|res| {abs(r2).max():.3f} MPa")
    # out-of-sample: fit on 500 and 200 only, predict the others
    k = [0, int(np.argmin(dx))]; c3 = np.polyfit(x[k], y[k], 1)
    print("predict from 500+200 only:", {int(a): round(float(np.polyval(c3, b)), 3) for a, b in zip(dx, x)}, "measured", {int(a): round(float(b), 3) for a, b in zip(dx, y)})

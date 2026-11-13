#! /usr/bin/env python3
"""Probe 1 figure. usage: fig_probe1.py <audit_dir> <out_png>
y-z section at x=0 (+/- one element) of the elements attached to the first 4 split pairs, for runs a60 (as built), a15 (as built), g60 (scratch guard).
Fill = which id of pairs 0-3 the element uses: slave (blue), master (orange); red edge = element on the footwall side that uses a master id (mis-wiring)."""
import sys, math
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
A, out = sys.argv[1], sys.argv[2]
runs = [('a60_dx500', 60, 'as built, 60 deg (dx=500 m)'), ('a15_dx500', 15, 'as built, 15 deg (test.tpv36 geometry, dx=500 m)'), ('g60_dx500', 60, 'scratch guard, 60 deg')]
fig, axs = plt.subplots(1, 3, figsize=(15, 5))
for ax, (nm, dip, title) in zip(axs, runs):
    t = math.tan(math.radians(dip)); c = math.cos(math.radians(dip))
    P = np.loadtxt(A + '/%s/audit_pairs.0000.txt' % nm); P = P[np.argsort(-P[:, 5])][:4]
    sl, ma = set(P[:, 1].astype(int)), set(P[:, 2].astype(int))
    E = np.loadtxt(A + '/%s/audit_elems.0000.txt' % nm)
    seen = set()
    for e in E:
        nodes = e[2:10].astype(int); yz = e[10:26].reshape(8, 2)
        us, um = bool(sl & set(nodes)), bool(ma & set(nodes))
        if not (us or um): continue
        pts = np.unique(np.round(yz, 3), axis=0)
        key = (tuple(map(tuple, pts)), int(e[1]), us, um)
        if key in seen: continue
        seen.add(key)
        d = (yz[:, 1] + yz[:, 0] * t) * c; off = d[np.abs(d) > 1.0]
        fw = off.size > 0 and (off < 0).all()
        mis = um and fw
        cx, cy = pts.mean(0)
        order = np.argsort(np.arctan2(pts[:, 1] - cy, pts[:, 0] - cx)) if len(pts) > 2 else range(len(pts))
        col = 'tab:orange' if um and not us else ('tab:blue' if us and not um else 'tab:purple')
        if len(pts) >= 3:
            ax.add_patch(Polygon(pts[list(order)] / 1e3, closed=True, fc=col, ec=('red' if mis else 'k'), lw=(2.5 if mis else 0.6), alpha=0.45))
        else:
            ax.plot(pts[:, 0] / 1e3, pts[:, 1] / 1e3, '-', color=col, lw=3)
    yy = np.linspace(0, 0.9 * 1e3 * (P[3, 4] / 1e3 + 0.2), 2)
    ax.plot([0, P[3, 4] * 1.15 / 1e3], [0, -P[3, 4] * 1.15 * t / 1e3], 'k--', lw=0.8)
    ax.plot(P[:, 4] / 1e3, P[:, 5] / 1e3, 'ko', ms=4)
    ax.set_aspect('equal'); ax.autoscale_view(); ax.set_xlabel('y (km)'); ax.set_ylabel('z (km)'); ax.set_title(title, fontsize=9)
    ax.autoscale(); ax.margins(0.1)
fig.suptitle('Elements attached to the trace node pair and the next three pairs (dots) at x=0; dashed = fault plane; blue slave-id, orange master-id, red edge = footwall element on a master id', fontsize=9)
fig.tight_layout(); fig.savefig(out, dpi=130)

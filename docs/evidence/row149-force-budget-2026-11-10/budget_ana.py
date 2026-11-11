#!/usr/bin/env python3
"""Row 149 force budget: reduce the diagEl/diagElS/diagNd dumps of the instrumented
EQdyna build to CSVs.

usage: budget_ana.py <runs_root> <morning_root> <out_dir>

Per split-node pair at x = 0 (rows 1..3 from the free surface) and per side:
  side traction of part p along e:  tM_p = +F_pM.e/arn,  tS_p = -F_pS.e/arn   [MPa]
  code metric (faulting.f90:93-95,145-155):
    T_e = (mS*mM*((vM-vS).e + [n only](uM-uS).e/dt)/dt + mS*FM.e - mM*FS.e)/((mS+mM)*arn) + T0*C_elastic
  so T_e = sum_p (mS*tM_p + mM*tS_p)/(mS+mM) + kinematic + T0*C_elastic.
Parts: 1 gravity, 2 mass damping, 3 B^T sigma, 4 stiffness damping, 5 hourglass (disp), 6 hourglass (rdampk*vel).
"""
import sys, os, csv, glob
import numpy as np

RUNS = [  # name, label, dx, C_elastic, kind, morning dir, morning dp000 n-stress at 4 s (MPa)
    ('t12_dx200', 'TPV12 60deg', 200, 1, 'dip', 't12_dx200_diag', -3.9865),
    ('t12_dx150', 'TPV12 60deg', 150, 1, 'dip', 't12_dx150_diag', -2.2977),
    ('t12v90_dx200', 'TPV12 dip=90 stand-in', 200, 1, 'vert', None, None),
    ('t12v90_dx150', 'TPV12 dip=90 stand-in', 150, 1, 'vert', None, None),
    ('t8_dx200', 'TPV8 vertical', 200, 1, 'vert', None, None),
    ('t8_dx150', 'TPV8 vertical', 150, 1, 'vert', None, None),
    ('t13_dx200', 'TPV13 60deg', 200, 0, 'dip', 't13_dx200_diag', 1.9555),
    ('t13_dx150', 'TPV13 60deg', 150, 0, 'dip', 't13_dx150_diag', 5.4215),
]
PARTS = ['grav', 'mdamp', 'BtS', 'kdamp', 'hg', 'hgdamp']
DIRS = ['n', 's', 'd']
T_EVAL = 4.0


def load(pattern):
    rows = []
    for fn in sorted(glob.glob(pattern)):
        if os.path.getsize(fn) == 0:
            continue
        rows.append(np.loadtxt(fn, ndmin=2))
    if not rows:
        raise SystemExit('no data for ' + pattern)
    return np.vstack(rows)


def station(path):
    """time and n-stress (MPa) columns of a SCEC on-fault station file."""
    t, s = [], []
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) < 8:
                continue
            try:
                v = [float(x) for x in p[:8]]
            except ValueError:
                continue
            t.append(v[0])
            s.append(v[7])
    if not t:
        raise SystemExit('no station rows in ' + path)
    return np.array(t), np.array(s)


def at_time(t, t0=T_EVAL):
    return int(np.argmin(np.abs(np.asarray(t) - t0)))


def reduce_run(rd):
    nd = load(os.path.join(rd, 'diagNd.*.txt'))
    el = load(os.path.join(rd, 'diagEl.*.txt'))
    es = load(os.path.join(rd, 'diagElS.*.txt'))
    # pairs: drop the second copy written by the neighbouring rank
    _, iu = np.unique(np.round(nd[:, [0, 4]], 3), axis=0, return_index=True)
    ncopies = len(nd) / len(iu)
    nd = nd[np.sort(iu)]
    # elements: dedup on (step, centroid, node xyz)
    _, ie = np.unique(np.round(el[:, [0, 4, 5, 6, 8, 9, 10]], 3), axis=0, return_index=True)
    el_dups = len(el) - len(ie)
    el = el[np.sort(ie)]
    _, ix = np.unique(np.round(es[:, [0, 4, 5, 6]], 3), axis=0, return_index=True)
    es = es[np.sort(ix)]
    elf = el[el[:, 11] > 0]
    F = {}
    for r in elf:
        k = (int(r[0]), int(r[11]), round(float(r[10]), 2))
        F[k] = F.get(k, 0.0) + r[13:31].reshape(6, 3)
    return dict(nd=nd, el=el, es=es, F=F, ncopies=ncopies, el_dups=el_dups)


def dot(a, b):
    return np.einsum('ij,ij->i', a, b)


def budget(red, cel):
    """Per pair row: part tractions per side and direction, kinematic term, closures."""
    nd, F = red['nd'], red['F']
    zs = np.flip(np.unique(np.round(nd[:, 4], 2)))            # row 1 = free surface
    rows = []
    for j, z in enumerate(zs):
        a = nd[np.abs(nd[:, 4] - z) < 0.01]
        a = a[np.argsort(a[:, 0])]
        n = len(a)
        arn, mS, mM, dt = a[:, 5], a[:, 8], a[:, 9], a[:, 49]
        e = dict(n=a[:, 10:13], s=a[:, 13:16], d=a[:, 16:19])
        FS = np.array([F[(int(nt), 1, round(float(z), 2))] for nt in a[:, 0]])   # (n, 6, 3)
        FM = np.array([F[(int(nt), 2, round(float(z), 2))] for nt in a[:, 0]])
        fpS, fpM, fqS, fqM = a[:, 19:22], a[:, 22:25], a[:, 25:28], a[:, 28:31]
        vS, vM = a[:, 31:34], a[:, 34:37]
        r = dict(row=j + 1, z=z, y=a[0, 3], t=a[:, 1], nt=a[:, 0], arn=arn, mS=mS, mM=mM)
        # closure 1: element parts summed over elements and ranks vs the assembled force
        scale = max(np.abs(fpS).max(), np.abs(fpM).max(), 1.0)
        r['clos_parts'] = max(np.abs(FS.sum(1) - fpS).max(), np.abs(FM.sum(1) - fpM).max()) / scale
        # closure 2: m*(v(n+1)-v(n))/dt vs the post-faulting force
        r['clos_inertia'] = float('nan')
        if n > 1:
            iS = mS[:-1, None] * (vS[1:] - vS[:-1]) / dt[:-1, None] - fqS[:-1]
            iM = mM[:-1, None] * (vM[1:] - vM[:-1]) / dt[:-1, None] - fqM[:-1]
            s2 = max(np.abs(fqS).max(), np.abs(fqM).max(), 1.0)
            r['clos_inertia'] = max(np.abs(iS).max(), np.abs(iM).max()) / s2
        budget_dirs(r, a, e, FS, FM, cel)
        rows.append(r)
    return rows


def budget_dirs(r, a, e, FS, FM, cel):
    arn, mS, mM, dt = a[:, 5], a[:, 8], a[:, 9], a[:, 49]
    fpS, fpM, fqS, fqM = a[:, 19:22], a[:, 22:25], a[:, 25:28], a[:, 28:31]
    vS, vM, uS, uM = a[:, 31:34], a[:, 34:37], a[:, 37:40], a[:, 40:43]
    w = mS + mM
    for di, d in enumerate(DIRS):
        ed = e[d]
        for p, pn in enumerate(PARTS):
            r['tM_' + pn + '_' + d] = dot(FM[:, p, :], ed) / arn / 1e6
            r['tS_' + pn + '_' + d] = -dot(FS[:, p, :], ed) / arn / 1e6
        dv = dot(vM - vS, ed)
        du = dot(uM - uS, ed)
        if d != 'n':
            du = 0.0 * du
        kin = mS * mM * (dv + du / dt) / dt / (w * arn)
        r['kin_' + d] = kin / 1e6
        t0 = a[:, 46 + di]
        r['T0_' + d] = t0 / 1e6
        r['Ttrial_' + d] = (kin + (mS * dot(fpM, ed) - mM * dot(fpS, ed)) / (w * arn) + t0 * cel) / 1e6
        r['Tcode_' + d] = a[:, 43 + di] / 1e6
        # traction applied by faulting on each side (post - pre force), same sign convention
        r['tMapp_' + d] = -dot(fqM - fpM, ed) / arn / 1e6
        r['tSapp_' + d] = dot(fqS - fpS, ed) / arn / 1e6
        # net (inertial) force per side as a traction: m*a = F_post
        r['tMnet_' + d] = dot(fqM, ed) / arn / 1e6
        r['tSnet_' + d] = -dot(fqS, ed) / arn / 1e6


def elem_stress(red, un, ud):
    """Top three element rows per side at each step: stress projections (MPa).
    side = +1 for the un side of the fault (hanging wall for the dipping cases)."""
    es = red['es']
    nd = red['nd']
    zf = nd[:, 4]
    yf = nd[:, 3]
    slope = 0.0
    if np.abs(zf).max() > 0:
        k = int(np.argmax(np.abs(zf)))
        slope = yf[k] / zf[k]
    zc = np.flip(np.unique(np.round(es[:, 6], 2)))
    out = []
    ez = np.array([0.0, 0.0, 1.0])
    for r in es:
        s = r[7:13]
        sig = np.array([[s[0], s[5], s[4]], [s[5], s[1], s[3]], [s[4], s[3], s[2]]])
        side = 1 if r[5] > slope * r[6] else -1
        row = int(np.argmin(np.abs(zc - r[6]))) + 1
        tn = sig @ un
        tz = sig @ ez
        out.append((r[0], r[1], row, side, r[4], tn @ un / 1e6, tn @ ud / 1e6, tz @ un / 1e6,
                    tz @ ud / 1e6, tz @ ez / 1e6, r[13], r[14], r[15]))
    return np.array(out)


def main():
    runs_root, morning_root, out_dir = sys.argv[1:4]
    os.makedirs(out_dir, exist_ok=True)
    ts_rows, at_rows, es_rows, cl_rows = [], [], [], []
    skeys = None
    for name, label, dx, cel, kind, mdir, mval in RUNS:
        rd = os.path.join(runs_root, name)
        if not os.path.exists(os.path.join(rd, 'DONE')):
            print('NOT RUN / not finished:', name)
            continue
        red = reduce_run(rd)
        rows = budget(red, cel)
        if skeys is None:
            skeys = [k for k in rows[0] if k[0] in 'tTk' and k not in ('t',)]
        st_t, st_s = station(os.path.join(rd, 'faultst000dp000.txt'))
        ist = at_time(st_t)
        for r in rows:
            i4 = at_time(r['t'])
            for i in range(len(r['t'])):
                if i % 2 == 0 or i == i4:
                    ts_rows.append([name, label, dx, r['row'], '%.3f' % r['z'], '%.6f' % r['t'][i]]
                                   + ['%.6g' % r[k][i] for k in skeys])
            at_rows.append([name, label, dx, r['row'], '%.3f' % r['z'], '%.6f' % r['t'][i4]]
                           + ['%.6g' % r[k][i4] for k in skeys])
            cl = [name, dx, r['row'], '%.3e' % r['clos_parts'], '%.3e' % r['clos_inertia'],
                  '%.3e' % np.abs(r['Ttrial_n'] - r['Tcode_n']).max(), '', '', '', '', '', '']
            if r['row'] == 1:
                cl[6:9] = ['%.6f' % st_t[ist], '%.4f' % st_s[ist], '%.4f' % r['Ttrial_n'][i4]]
                cl[9] = len(st_t)
                if mdir:
                    mt, ms = station(os.path.join(morning_root, mdir, 'faultst000dp000.txt'))
                    cl[10:12] = ['%.4f' % ms[at_time(mt)], len(mt)]
            cl_rows.append(cl)
        a0 = red['nd'][0]
        es = elem_stress(red, a0[10:13], a0[16:19])
        t4 = es[at_time(es[:, 1]), 1]
        for q in es[(es[:, 1] == t4) | (es[:, 0] == es[:, 0].min())]:
            es_rows.append([name, label, dx] + ['%.6g' % v for v in q])
        print('%-14s pairs-copies %.1f elem-dups %d  steps %d' % (name, red['ncopies'], red['el_dups'], len(rows[0]['t'])))
    head = ['run', 'label', 'dx_m', 'row', 'z_m', 't_s']
    write_csv(os.path.join(out_dir, 'timeseries.csv'), head + skeys, ts_rows)
    write_csv(os.path.join(out_dir, 'at4s.csv'), head + skeys, at_rows)
    write_csv(os.path.join(out_dir, 'elem_stress.csv'),
              ['run', 'label', 'dx_m', 'step', 't_s', 'row', 'side', 'x_centroid_m', 'sig_nn', 'sig_nd', 'sig_zn',
               'sig_zd', 'sig_zz', 'rjust', 'pstrain', 'vol_m3'], es_rows)
    write_csv(os.path.join(out_dir, 'closure.csv'),
              ['run', 'dx_m', 'row', 'clos_parts_rel', 'clos_inertia_rel', 'max_abs_Ttrial_minus_Tcode_MPa',
               't_station_s', 'station_nstress_MPa', 'Ttrial_n_MPa', 'station_rows', 'morning_station_nstress_MPa',
               'morning_station_rows'], cl_rows)


def write_csv(path, head, rows):
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(head)
        w.writerows(rows)
    with open(path, newline='') as f:
        back = list(csv.reader(f))
    if len(back) != len(rows) + 1 or any(len(b) != len(head) for b in back):
        raise SystemExit('csv re-read mismatch: ' + path)
    print('wrote %s: %d rows x %d cols (re-read ok)' % (path, len(rows), len(head)))


if __name__ == '__main__':
    main()

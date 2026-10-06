#! /usr/bin/env python3
"""
Compares a test.tpv34 run with EQdyna's own 2016 SCEC TPV34 submission
(50 m, 20 s, EQdyna3d v3.2.3, author D.Liu), FROM THE ARCHIVED DATA -- the
rule 17 step 6 independent check for test.tpv34, committed as a script so
its numbers carry provenance instead of living in prose (PROJECT_RULES rule
4 / rule 6; the TPV29 lesson, pathway item 28).

REPORT-ONLY, same as evidence_tpv29/35_scec_comparison.py: never a gate,
never wired into testsys/run.py, invoked by hand. Prints the numbers, says
what each one tests and what it cannot test, writes a JSON snapshot, and
always exits 0.

TWO KINDS OF RUN DIRECTORY
--------------------------
The run's own dz and term are read off a station-file header, never
assumed, so the script takes either:

  * the GATE configuration (500 m, 5 s, 4 ranks, the shipped 500 m grid) --
    what the case is gated at; only the early-time items below mean much;
  * the VALIDATION run the owner accepted on 2026-10-06 (board row 19(c),
    PR #97: "100 m is good enough for validation. Show me the overlay of
    cplot and ts"): dx = 100 m, term = 20 s (the spec's full term), on a
    100 m CVM-H element-centre grid RE-EXTRACTED with the committed
    case_input/test.tpv34/extract_cvmh_grid.py recipe (a scratch case: the
    committed case only admits its shipped 500 m grid, and the gate
    reference is untouched). The comparison is then 100 m against the 50 m
    submission -- NOT matched resolution; the owner accepted that.

Recipe for the validation run (scratch, never committed: the 100 m grid is
49.9 M rows): extract the 100 m grid with extract_cvmh_grid.py's
cell_centres / to_utm / run_vx_lite / Step 6-7 logic (chunked, stored as
.npy); create.newcase <dir> test.tpv34; in <dir> point tpv34Tools'
GRID_FILE/SHIPPED_DX at it and set par.dx = 100, par.term = 20 and the rank
layout; python3 case.setup; mpirun -np N eqdyna.

WHAT IS COMPARED
----------------
The 2016 run used CVM-H 15.1.0, this repo's extractions use 15.1.1 (see
case_input/test.tpv34/README.md).

  1. INITIAL STRESS at the 35 on-fault stations: first-sample h-shear-stress
     and n-stress, run vs 2016 (the 2016 files start at t = 0.004 s). Tests
     the whole Part 2 + Part 3 chain -- CVM-H extraction, UTM mapping,
     clamp, mu/mu0 scaling, nucleation increment -- against an independent
     2016 extraction. A sign flip or tens-of-percent miss would mean a wrong
     frame, UTM or scaling.
  2. RUPTURE ARRIVAL at the on-fault stations: first time h-slip-rate
     exceeds SCEC's 1 mm/s threshold, t_run - t_2016.
  3. RUPTURE-TIME FIELD over the whole fault: the run's frt (deduped) vs the
     2016 cplot AT THE NODES THEY SHARE -- the 2016 grid is 50 m, so every
     node of a 100 m (or 500 m) run is a 2016 node: no interpolation.
     Reports the fraction ruptured by the term on each run's own grid, the
     fraction of shared nodes ruptured in both / in one only, and the
     signed median, median |dt| and p90 |dt| of t_run - t_2016 where both
     ruptured.
  4. PEAK SLIP RATE and FINAL SLIP at the on-fault stations both runs
     ruptured within the term: run/2016 ratio of max h-slip-rate within the
     term and of h-slip at the term.
  5. OFF-FAULT PEAK VELOCITY at the 56 off-fault stations: run/2016 ratio
     of peak |h-vel| and |v-vel| within the term (ratio stats only where the
     2016 peak is >= 1% of the largest 2016 peak, as the TPV35 script does).
  6. --plot-dir: the owner's overlays. cplot: rupture-time contours on the
     fault, run vs 2016, same axes, same levels. ts: on-fault (h-slip,
     h-slip-rate, h-shear-stress) and off-fault (h-, v-, n-vel) series at a
     representative station set (ON_FAULT_TS / OFF_FAULT_TS below).

The 2016 archive's 2.4 km stations are named dp024; at 100 m the run's are
too, at 500 m they are dp025 (paired by the run's rounding rule, printed).

Writes a timestamped JSON snapshot to testsys/parity/evidence_output/
(gitignored), like the TPV29/TPV35 scripts.
"""
import argparse
import datetime
import glob
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TESTSYS = os.path.dirname(HERE)
REPO_ROOT = os.path.dirname(TESTSYS)
OUT_DIR = os.path.join(HERE, 'evidence_output')
SUBMISSION_DIR_DEFAULT = os.path.join(REPO_ROOT, 'scec_archive', 'tpv34',
                                      'eqdyna3d-v3.2.3-50m-2016')
SLIP_RATE_THRESHOLD = 1.0e-3   # m/s, SCEC rupture-front definition
UNRUPTURED_CPLOT = 1000.0      # s, the 2016 cplot sentinel
UNRUPTURED_FRT = 99999.0       # s, EQdyna frt column 4 sentinel
FAULT_AREA_KM2 = 30.0 * 15.0
ON_FAULT_X_KM = (-12.0, -6.0, 0.0, 6.0, 12.0)
ON_FAULT_DEPTH_KM_SPEC = (0.0, 1.0, 2.4, 5.0, 7.5, 10.0, 12.0)
OFF_FAULT_RATIO_FLOOR = 0.01
# Representative ts sets: (x_km, spec depth_km) on the fault -- hypocentre,
# both rupture directions at the surface, mid-depth and deep -- and off-fault
# (normal_km, x_km, depth_km): near/far field, both sides, both depths, and
# the in-plane station beyond the fault end.
ON_FAULT_TS = ((0.0, 7.5), (0.0, 0.0), (-6.0, 2.4), (6.0, 2.4), (-12.0, 5.0),
               (12.0, 5.0), (-6.0, 10.0), (6.0, 12.0))
OFF_FAULT_TS = ((3.0, 0.0, 0.0), (-3.0, 10.0, 0.0), (9.0, -10.0, 0.0),
                (-9.0, 20.0, 0.0), (15.0, 0.0, 0.0), (0.0, 20.0, 0.0),
                (3.0, -10.0, 2.4), (-15.0, -15.0, 2.4))


def _numeric_rows(path):
    rows = []
    with open(path) as f:
        for line in f:
            s = line.split()
            if not s or s[0].startswith('#'):
                continue
            try:
                rows.append([float(v) for v in s])
            except ValueError:
                continue
    return np.array(rows)


def _hm(km):
    """Station-name field in hectometres, as both writers format it."""
    v = int(round(km * 10.0))
    return '%03d' % v if v >= 0 else '-%03d' % (-v)


def _depth_for_run(d_km, run_dz_m):
    return round(d_km * 1.0e3 / run_dz_m) * run_dz_m / 1.0e3


def station_pairs(run_dz_m):
    """[(run_name, archive_name, x_km, depth_km_spec, depth_km_run), ...]"""
    out = []
    for x in ON_FAULT_X_KM:
        for d in ON_FAULT_DEPTH_KM_SPEC:
            dg = _depth_for_run(d, run_dz_m)
            out.append(('faultst%sdp%s.txt' % (_hm(x), _hm(dg)),
                        'faultst%sdp%s' % (_hm(x), _hm(d)), x, d, dg))
    return out


def off_fault_names():
    """[(run_name, archive_name, normal_km, x_km, depth_km)] for the spec's
    56 off-fault stations (tpv34Tools.offFaultStationsKm order). The run
    names them from the REQUESTED coordinate (library_output.f90), so the
    names are resolution-independent and equal the archive's plus .txt."""
    out = []
    for depth in (0.0, 2.4):
        pts = [(z, x) for z in (-9.0, -3.0, 3.0, 9.0) for x in (-20.0, -10.0, 0.0, 10.0, 20.0)]
        pts += [(z, x) for z in (-15.0, 15.0) for x in (-15.0, 0.0, 15.0)]
        pts += [(0.0, x) for x in (-20.0, 20.0)]
        for z, x in pts:
            name = 'body%sst%sdp%s' % (_hm(z), _hm(x), _hm(depth))
            out.append((name + '.txt', name, z, x, depth))
    return out


def first_sample_stress(a):
    return float(a[0, 3]), float(a[0, 7]), float(a[0, 0])


def arrival_time(a, term):
    a = a[a[:, 0] <= term + 1e-9]
    hit = np.nonzero(a[:, 2] > SLIP_RATE_THRESHOLD)[0]
    return float(a[hit[0], 0]) if hit.size else None


def at_time(a, t, col):
    i = int(np.argmin(np.abs(a[:, 0] - t)))
    return float(a[i, col]), float(a[i, 0])


def peak_abs(a, term, col):
    a = a[a[:, 0] <= term + 1e-9]
    i = int(np.argmax(np.abs(a[:, col])))
    return float(abs(a[i, col])), float(a[i, 0])


def read_cplot(path):
    """(x_m, depth_m, t_s) of the 2016 cplot, sentinel kept."""
    a = _numeric_rows(path)
    return a[:, 0], a[:, 1], a[:, 2]


def read_frt(run_dir):
    """(x_m, depth_m, t_s) of the run's frt files, deduped on rounded xyz."""
    files = sorted(glob.glob(os.path.join(run_dir, 'frt.txt*')))
    if not files:
        raise SystemExit('%s has no frt.txt* files' % run_dir)
    a = np.vstack([np.loadtxt(f, usecols=(0, 1, 2, 3)) for f in files])
    _, first = np.unique(np.round(a[:, :3], 3), axis=0, return_index=True)
    a = a[first]
    return a[:, 0], -a[:, 2], a[:, 3]


def run_term_and_dz(run_dir):
    """Read the run's term (s) and dz (m) off an on-fault station header."""
    path = os.path.join(run_dir, 'faultst000dp075.txt')
    dz = term = dt = None
    with open(path) as f:
        for line in f:
            key = line.lstrip('# ').split('=')[0].strip()
            if key == 'element_size':
                dz = float(line.split('=')[1])
            elif key == 'time_step':       # 'num_time_steps' also contains this substring
                dt = float(line.split('=')[1].split()[0])
            elif key == 'num_time_steps':
                term = dt * int(line.split('=')[1])
    return term, dz


def provenance():
    def git(*args):
        try:
            return subprocess.check_output(['git'] + list(args), cwd=REPO_ROOT,
                                           text=True).strip()
        except Exception as ex:      # noqa: BLE001 -- provenance is best effort
            return 'unavailable (%s)' % ex
    return dict(git_sha=git('rev-parse', 'HEAD'), git_dirty=bool(git('status', '--porcelain')),
                when_utc=datetime.datetime.utcnow().isoformat(timespec='seconds'),
                host=os.uname().nodename)


def _stats(v):
    v = np.asarray(v, float)
    if not v.size:
        return dict(n=0)
    return dict(n=int(v.size), median=float(np.median(v)), min=float(v.min()),
                max=float(v.max()), p10=float(np.percentile(v, 10)),
                p90=float(np.percentile(v, 90)))


def rupture_field(sub, run, term):
    """Item 3: whole-fault rupture time at the shared nodes."""
    cx, cd, ct = read_cplot(os.path.join(sub, 'cplot'))
    fx, fd, ft = read_frt(run)
    c_rupt = (ct < UNRUPTURED_CPLOT) & (ct <= term)
    f_rupt = (ft < UNRUPTURED_FRT) & (ft <= term)
    key = lambda x, d: np.round(x).astype(np.int64) * 100000 + np.round(d).astype(np.int64)
    ck, fk = key(cx, cd), key(fx, fd)
    _, ic, jf = np.intersect1d(ck, fk, return_indices=True)
    both = c_rupt[ic] & f_rupt[jf]
    dt = ft[jf][both] - ct[ic][both]
    return dict(
        ruptured_fraction_2016=float(c_rupt.mean()), nodes_2016=int(ct.size),
        ruptured_fraction_run=float(f_rupt.mean()), nodes_run=int(ft.size),
        shared_nodes=int(ic.size),
        shared_ruptured_both=int(both.sum()),
        shared_ruptured_run_only=int((f_rupt[jf] & ~c_rupt[ic]).sum()),
        shared_ruptured_2016_only=int((c_rupt[ic] & ~f_rupt[jf]).sum()),
        dt_signed_median_s=float(np.median(dt)) if dt.size else None,
        dt_abs_median_s=float(np.median(np.abs(dt))) if dt.size else None,
        dt_abs_p90_s=float(np.percentile(np.abs(dt), 90)) if dt.size else None,
        dt_abs_max_s=float(np.abs(dt).max()) if dt.size else None)


def compare(sub, run):
    term, dz = run_term_and_dz(run)
    r = dict(run_term_s=term, run_dz_m=dz, stations=[], off_fault=[])
    for rname, aname, x, d, dg in station_pairs(dz):
        rp, apath = os.path.join(run, rname), os.path.join(sub, aname)
        if not (os.path.isfile(rp) and os.path.isfile(apath)):
            r['stations'].append(dict(run=rname, archive=aname, missing=True))
            continue
        b, a = _numeric_rows(rp), _numeric_rows(apath)
        gs, gn, gt0 = first_sample_stress(b)
        as_, an, at0 = first_sample_stress(a)
        psr_r, tpsr_r = peak_abs(b, term, 2)
        psr_a, tpsr_a = peak_abs(a, term, 2)
        sl_r, _ = at_time(b, term, 1)
        sl_a, _ = at_time(a, term, 1)
        r['stations'].append(dict(
            run=rname, archive=aname, x_km=x, depth_km_spec=d, depth_km_run=dg,
            shear0_run_MPa=gs, shear0_2016_MPa=as_, shear0_ratio=gs / as_,
            normal0_run_MPa=gn, normal0_2016_MPa=an, normal0_ratio=gn / an,
            first_sample_s=(gt0, at0),
            arrival_run_s=arrival_time(b, term), arrival_2016_s=arrival_time(a, term),
            peak_sr_run=psr_r, peak_sr_2016=psr_a,
            peak_sr_ratio=(psr_r / psr_a if psr_a > SLIP_RATE_THRESHOLD else None),
            slip_term_run_m=sl_r, slip_term_2016_m=sl_a,
            slip_ratio=(sl_r / sl_a if abs(sl_a) > 1e-3 else None)))
    ok = [s for s in r['stations'] if not s.get('missing')]
    sr = np.array([s['shear0_ratio'] for s in ok])
    nr = np.array([s['normal0_ratio'] for s in ok])
    r['initial_stress'] = dict(
        n=len(ok), shear_ratio=_stats(sr), normal_ratio=_stats(nr),
        sign_agreement=bool(np.all(sr > 0) and np.all(nr > 0)))
    both = [s for s in ok if s['arrival_run_s'] is not None and s['arrival_2016_s'] is not None]
    dts = np.array([s['arrival_run_s'] - s['arrival_2016_s'] for s in both])
    r['arrival'] = dict(
        n_both=len(both), n_run_only=sum(1 for s in ok if s['arrival_run_s'] is not None
                                         and s['arrival_2016_s'] is None),
        n_2016_only=sum(1 for s in ok if s['arrival_2016_s'] is not None
                        and s['arrival_run_s'] is None),
        dt_signed=_stats(dts), dt_abs=_stats(np.abs(dts)))
    r['rupture_field'] = rupture_field(sub, run, term)
    # Ratios only where BOTH runs ruptured the station within the term: a
    # station only one run reached is counted in item 2, not averaged here as 0.
    r['peak_slip_rate_ratio'] = _stats([s['peak_sr_ratio'] for s in both if s['peak_sr_ratio']])
    r['final_slip_ratio'] = _stats([s['slip_ratio'] for s in both if s['slip_ratio']])
    hyp = [s for s in ok if s['archive'] == 'faultst000dp075']
    r['hypocentre'] = hyp[0] if hyp else None

    for rname, aname, z, x, d in off_fault_names():
        rp, apath = os.path.join(run, rname), os.path.join(sub, aname)
        if not (os.path.isfile(rp) and os.path.isfile(apath)):
            r['off_fault'].append(dict(run=rname, archive=aname, missing=True,
                                       run_present=os.path.isfile(rp),
                                       archive_present=os.path.isfile(apath)))
            continue
        b, a = _numeric_rows(rp), _numeric_rows(apath)
        e = dict(run=rname, archive=aname, normal_km=z, x_km=x, depth_km=d)
        for col, q in ((2, 'hvel'), (4, 'vvel')):
            pr, tr = peak_abs(b, term, col)
            pa, ta = peak_abs(a, term, col)
            e.update({'peak_%s_run' % q: pr, 'peak_%s_2016' % q: pa,
                      'peak_%s_ratio' % q: (pr / pa if pa else None),
                      'peak_%s_dt_s' % q: tr - ta})
        r['off_fault'].append(e)
    ofok = [s for s in r['off_fault'] if not s.get('missing')]
    summ = dict(n=len(ofok), n_missing=len(r['off_fault']) - len(ofok),
                ratio_floor_fraction_of_max=OFF_FAULT_RATIO_FLOOR)
    for q in ('hvel', 'vvel'):
        pmax = max([s['peak_%s_2016' % q] for s in ofok] or [0.0])
        use = [s for s in ofok if s['peak_%s_2016' % q] >= OFF_FAULT_RATIO_FLOOR * pmax > 0]
        summ['%s_ratio' % q] = _stats([s['peak_%s_ratio' % q] for s in use])
        summ['%s_peak_time_abs_dt_s' % q] = _stats([abs(s['peak_%s_dt_s' % q]) for s in use])
    r['off_fault_summary'] = summ
    return r


def _fmt_stats(s, f='%.3f'):
    if not s.get('n'):
        return 'n=0'
    return ('n=%d median ' + f + '  p10 ' + f + '  p90 ' + f + '  range ' + f + ' .. ' + f) % (
        s['n'], s['median'], s['p10'], s['p90'], s['min'], s['max'])


def print_report(r, sub, run):
    print('TPV34: run %s (dz %g m, term %g s) vs 2016 SCEC submission %s (50 m, 20 s)'
          % (run, r['run_dz_m'], r['run_term_s'], sub))
    print('NOT matched resolution: %g m run vs the 50 m submission (owner-accepted 2026-10-06, '
          'row 19(c)); CVM-H 15.1.1 vs 15.1.0.' % r['run_dz_m'])
    i = r['initial_stress']
    print('\n1. initial stress, run/2016 ratio over %d on-fault stations (signs agree: %s):'
          % (i['n'], i['sign_agreement']))
    print('   h-shear  %s' % _fmt_stats(i['shear_ratio'], '%.4f'))
    print('   n-stress %s' % _fmt_stats(i['normal_ratio'], '%.4f'))
    print('   %-22s %-18s %8s %8s %8s %8s %7s %7s %7s %7s %6s %6s' % (
        'run', '2016', 'tau_r', 'tau_16', 'sn_r', 'sn_16', 't_r', 't_16', 'psr_r', 'psr_16',
        'D_r', 'D_16'))
    nn = lambda v: float('nan') if v is None else v
    for s in r['stations']:
        if s.get('missing'):
            print('   %-22s %-18s MISSING' % (s['run'], s['archive']))
            continue
        print('   %-22s %-18s %8.3f %8.3f %8.3f %8.3f %7.3f %7.3f %7.3f %7.3f %6.3f %6.3f' % (
            s['run'], s['archive'], s['shear0_run_MPa'], s['shear0_2016_MPa'],
            s['normal0_run_MPa'], s['normal0_2016_MPa'], nn(s['arrival_run_s']),
            nn(s['arrival_2016_s']), s['peak_sr_run'], s['peak_sr_2016'],
            s['slip_term_run_m'], s['slip_term_2016_m']))
    a = r['arrival']
    print('\n2. station rupture arrival (h-slip-rate > 1 mm/s) within %g s: %d in both, '
          '%d run-only, %d 2016-only' % (r['run_term_s'], a['n_both'], a['n_run_only'],
                                         a['n_2016_only']))
    print('   run - 2016 signed: %s' % _fmt_stats(a['dt_signed']))
    print('   |run - 2016|     : %s' % _fmt_stats(a['dt_abs']))
    f = r['rupture_field']
    print('\n3. rupture-time field by %g s: ruptured fraction 2016 %.4f (%d nodes, 50 m), '
          'run %.4f (%d nodes)' % (r['run_term_s'], f['ruptured_fraction_2016'],
                                   f['nodes_2016'], f['ruptured_fraction_run'], f['nodes_run']))
    print('   shared nodes %d: ruptured in both %d, run only %d, 2016 only %d'
          % (f['shared_nodes'], f['shared_ruptured_both'], f['shared_ruptured_run_only'],
             f['shared_ruptured_2016_only']))
    if f['dt_abs_median_s'] is not None:
        print('   t_run - t_2016 where both ruptured: signed median %+.3f s, median |dt| %.3f s, '
              'p90 |dt| %.3f s, max |dt| %.3f s' % (f['dt_signed_median_s'], f['dt_abs_median_s'],
                                                    f['dt_abs_p90_s'], f['dt_abs_max_s']))
    print('\n4. on-fault stations ruptured in both runs, run/2016:')
    print('   peak h-slip-rate ratio  %s' % _fmt_stats(r['peak_slip_rate_ratio']))
    print('   final h-slip ratio      %s' % _fmt_stats(r['final_slip_ratio']))
    o = r['off_fault_summary']
    print('\n5. off-fault peak velocity, run/2016 over %d of %d stations (%d missing), '
          'stats where 2016 peak >= %g x max:' % (o['n'], o['n'] + o['n_missing'], o['n_missing'],
                                                  o['ratio_floor_fraction_of_max']))
    print('   peak |h-vel| ratio %s' % _fmt_stats(o['hvel_ratio']))
    print('   peak |v-vel| ratio %s' % _fmt_stats(o['vvel_ratio']))
    print('   |t_peak dt| h-vel  %s' % _fmt_stats(o['hvel_peak_time_abs_dt_s']))
    for s in r['off_fault']:
        if s.get('missing'):
            print('   %-22s MISSING (run %s, archive %s)' % (s['run'], s['run_present'],
                                                           s['archive_present']))
    print('\nThis is a report (%g m / %g s vs 50 m / 20 s, CVM-H 15.1.1 vs 15.1.0), not a gate.'
          % (r['run_dz_m'], r['run_term_s']))


def make_overlays(r, sub, run, plot_dir):
    """Owner 2026-10-06 (board row 19(c)): cplot and ts overlays, run vs 2016.
    Writes three PNGs into plot_dir and returns their absolute paths."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    os.makedirs(plot_dir, exist_ok=True)
    term, dz = r['run_term_s'], r['run_dz_m']
    out = []

    # (1) cplot: rupture-time contours, same axes, same levels.
    cx, cd, ct = read_cplot(os.path.join(sub, 'cplot'))
    xs, ks = np.unique(cx), np.unique(cd)
    T = np.full((ks.size, xs.size), np.nan)
    T[np.searchsorted(ks, cd), np.searchsorted(xs, cx)] = ct
    T[(T >= UNRUPTURED_CPLOT) | (T > term)] = np.nan
    fx, fd, ft = read_frt(run)
    m = (ft < UNRUPTURED_FRT) & (ft <= term)
    levels = np.arange(1.0, term + 1e-9, 1.0)
    fig, ax = plt.subplots(figsize=(11, 6))
    ax.contour(xs / 1e3, ks / 1e3, T, levels=levels, colors='k', linewidths=0.9)
    ax.tricontour(fx[m] / 1e3, fd[m] / 1e3, ft[m], levels=levels, colors='r',
                  linewidths=0.9, linestyles='--')
    for x in ON_FAULT_X_KM:
        for d in ON_FAULT_DEPTH_KM_SPEC:
            ax.plot(x, d, 'b^', ms=4)
    ax.plot(0, 7.5, 'k*', ms=12)
    ax.set_xlim(-15, 15); ax.set_ylim(15, 0)
    ax.set_aspect('equal')
    ax.set_xlabel('along strike (km)')
    ax.set_ylabel('down dip (km)')
    f = r['rupture_field']
    ax.set_title('TPV34 rupture-time contours every 1 s to %g s: 2016 submission 50 m (black) vs '
                 'this run %g m (red dashed)\nruptured fraction 2016 %.3f, run %.3f; '
                 '|dt| median %.3f s, p90 %.3f s (NOT matched resolution)'
                 % (term, dz, f['ruptured_fraction_2016'], f['ruptured_fraction_run'],
                    f['dt_abs_median_s'] or float('nan'), f['dt_abs_p90_s'] or float('nan')),
                 fontsize=9)
    ax.legend([Line2D([], [], color='k'), Line2D([], [], color='r', ls='--'),
               Line2D([], [], color='b', marker='^', ls='')],
              ['2016 submission (50 m)', 'this run (%g m)' % dz, 'on-fault stations'],
              fontsize=8, loc='lower left')
    pth = os.path.abspath(os.path.join(plot_dir, 'tpv34_cplot_overlay.png'))
    fig.tight_layout(); fig.savefig(pth, dpi=150); plt.close(fig); out.append(pth)

    # (2) on-fault ts.
    cols = ((1, 'h-slip (m)'), (2, 'h-slip-rate (m/s)'), (3, 'h-shear-stress (MPa)'))
    fig, axs = plt.subplots(len(ON_FAULT_TS), 3, figsize=(14, 1.9 * len(ON_FAULT_TS)), sharex=True)
    for i, (x, d) in enumerate(ON_FAULT_TS):
        an = 'faultst%sdp%s' % (_hm(x), _hm(d))
        rn = 'faultst%sdp%s.txt' % (_hm(x), _hm(_depth_for_run(d, dz)))
        a, b = _numeric_rows(os.path.join(sub, an)), _numeric_rows(os.path.join(run, rn))
        for j, (col, lab) in enumerate(cols):
            axs[i, j].plot(a[:, 0], a[:, col], 'k', lw=1, label='2016 50 m')
            axs[i, j].plot(b[:, 0], b[:, col], 'r--', lw=1, label='run %g m' % dz)
            axs[i, j].set_ylabel(lab, fontsize=8)
            axs[i, j].set_title('%s (x=%g km, depth %g km)' % (an, x, d), fontsize=8)
            axs[i, j].tick_params(labelsize=7)
    axs[0, 0].legend(fontsize=7)
    for a_ in axs[-1]:
        a_.set_xlabel('t (s)'); a_.set_xlim(0, term)
    fig.tight_layout()
    pth = os.path.abspath(os.path.join(plot_dir, 'tpv34_ts_onfault_overlay.png'))
    fig.savefig(pth, dpi=110); plt.close(fig); out.append(pth)

    # (3) off-fault ts.
    cols = ((2, 'h-vel (m/s)'), (4, 'v-vel (m/s)'), (6, 'n-vel (m/s)'))
    fig, axs = plt.subplots(len(OFF_FAULT_TS), 3, figsize=(14, 1.9 * len(OFF_FAULT_TS)), sharex=True)
    for i, (z, x, d) in enumerate(OFF_FAULT_TS):
        an = 'body%sst%sdp%s' % (_hm(z), _hm(x), _hm(d))
        a, b = _numeric_rows(os.path.join(sub, an)), _numeric_rows(os.path.join(run, an + '.txt'))
        for j, (col, lab) in enumerate(cols):
            axs[i, j].plot(a[:, 0], a[:, col], 'k', lw=1, label='2016 50 m')
            axs[i, j].plot(b[:, 0], b[:, col], 'r--', lw=1, label='run %g m' % dz)
            axs[i, j].set_ylabel(lab, fontsize=8)
            axs[i, j].set_title('%s (normal %g km, x %g km, depth %g km)' % (an, z, x, d), fontsize=8)
            axs[i, j].tick_params(labelsize=7)
    axs[0, 0].legend(fontsize=7)
    for a_ in axs[-1]:
        a_.set_xlabel('t (s)'); a_.set_xlim(0, term)
    fig.tight_layout()
    pth = os.path.abspath(os.path.join(plot_dir, 'tpv34_ts_offfault_overlay.png'))
    fig.savefig(pth, dpi=110); plt.close(fig); out.append(pth)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--submission-dir', default=SUBMISSION_DIR_DEFAULT)
    ap.add_argument('--run-dir', required=True,
                    help='a COMPLETED test.tpv34 run directory (faultst*.txt, body*.txt, '
                         'frt.txt*): the gate config, or the 100 m / 20 s validation run '
                         '(see the module docstring for its recipe)')
    ap.add_argument('--plot-dir', default=None,
                    help='write the cplot and ts overlay PNGs here (owner 2026-10-06, row 19(c))')
    ap.add_argument('--no-json', action='store_true')
    args = ap.parse_args()
    if not os.path.isdir(args.submission_dir):
        raise SystemExit('%s missing -- fetch per scec_archive/tpv34/*/PROVENANCE.md'
                         % args.submission_dir)
    r = compare(args.submission_dir, args.run_dir)
    print_report(r, args.submission_dir, args.run_dir)
    if args.plot_dir:
        r['overlay_pngs'] = make_overlays(r, args.submission_dir, args.run_dir, args.plot_dir)
        for pth in r['overlay_pngs']:
            print('overlay: %s' % pth)
    if not args.no_json:
        os.makedirs(OUT_DIR, exist_ok=True)
        stamp = datetime.datetime.utcnow().strftime('%Y%m%dT%H%M%SZ')
        path = os.path.join(OUT_DIR, 'tpv34_scec_comparison_%s.json' % stamp)
        with open(path, 'w') as f:
            json.dump(dict(provenance=provenance(), submission_dir=args.submission_dir,
                           run_dir=args.run_dir, result=r), f, indent=1)
        print('snapshot: %s' % path)
    return 0


if __name__ == '__main__':
    sys.exit(main())

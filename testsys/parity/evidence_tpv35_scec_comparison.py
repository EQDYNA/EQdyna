#! /usr/bin/env python3
"""
Compares a test.tpv35 run at the SPEC resolution/term (100 m, 18 s) with
EQdyna's own 2017 SCEC TPV35 submission (also 100 m, 18 s, same code lineage
-- author D.Liu, EQdyna3D v4.1.3) FROM THE ARCHIVED DATA -- the rule 17 step
6 independent check for test.tpv35, committed as a script so its numbers
carry provenance instead of living in prose (PROJECT_RULES rule 4 / rule 6;
the TPV29 lesson, pathway item 28). Modeled closely on
evidence_tpv34_scec_comparison.py's structure and conventions.

REPORT-ONLY: never a gate, never wired into testsys/run.py, invoked by hand.
Prints the numbers, says what each one tests and what it cannot test, writes
a JSON snapshot, and always exits 0.

WHAT IS DIFFERENT FROM TPV34's SCRIPT
--------------------------------------
TPV35's 2017 archive is at the SAME resolution (100 m) and SAME term (18 s)
as the spec recommends -- unlike TPV34's 2016 archive (50 m, 20 s) against a
500 m/5 s gate run. This script is meant to be run against a run-dir
produced at dx=100, term=18 (a SCRATCH case, NOT the committed 500 m/5 s
gate config -- see case_input/test.tpv35/user_defined_params.py's own
comment: "par.dx below is the gate value"). It still reads the run's actual
term/dz from a station-file header rather than hardcoding (same pattern as
TPV34's script), so it also tolerates being pointed at the gate run, in
which case every comparison below is reported with the mismatch plainly
visible in gate_term_s/gate_dz_m rather than hidden.

TPV35's off-fault stations are named (e.g. "4064_donna"), not numbered, and
the archive's on-fault files already use this project's own hectometre
station-naming convention (`faultst<x_hm>dp081`) -- the archive was produced
by this same codebase in 2017, so no cross-code renaming/pairing logic is
needed for on-fault stations; off-fault stations are paired by NAME, which
the archive shares with tpv35_station_locations.txt, against the run's own
body<y_hm>st<x_hm>dp000.txt files (library_output.f90:217-221), which are
not name-stamped and have to be re-derived from the station list's coordinates.

WHAT IS COMPARED
----------------
  1. INITIAL STRESS at the 7 on-fault stations (x = -25..+5 km step 5, depth
     8.1 km): first-sample h-shear-stress and n-stress, run vs 2017 archive.
     Tests the whole chain this case is new on: the official 100 m
     mu_s/tau0 grid decimation (tpv35Tools.faultGridForCase), the
     near/far-side 1D material assignment, and sign conventions -- against
     an independent 2017 run of the SAME data by the SAME code lineage.
     Expected near-exact agreement (same grid, same physics, same code
     family) if nothing has drifted; a tens-of-percent miss or sign flip
     means a real divergence, not resolution noise (unlike TPV34/TPV29).
  2. RUPTURE ARRIVAL at the on-fault stations that rupture within the run's
     term in BOTH runs: first time h-slip-rate exceeds SCEC's 1 mm/s
     threshold, t_run - t_2017.
  3. RUPTURED AREA AT THE TERM: fraction of the 40 x 15.5 km fault with
     rupture time <= term from the 2017 cplot (100 m nodes) vs the run's
     frt files, both node-counted on their own grids -- same grid if the
     run is also at 100 m.
  4. HYPOCENTRE SLIP AT THE TERM: h-slip at t = term at faultst000dp081
     (x = 0, the hypocentre's along-strike position).
  5. OFF-FAULT GROUND MOTION at the 42 named surface stations the archive
     ships (of 43 in the official list -- one, 4134_vyc, drops at this
     mesh; case.setup's own stderr already reports any station a run drops
     for the SAME reason): peak |h-vel| and peak |v-vel|, run vs 2017,
     ratio and arrival-time difference of the peak. This is the ground-
     motion half of a VALIDATION benchmark (TPV35 is meant to be compared
     to real Parkfield recordings, not just another code) -- comparing
     against a real-code archive at matched resolution is the nearest
     thing to that this script can do; it does NOT compare to the actual
     NGA-West2 seismic recordings Ma et al. (2008) fit to, which is a
     separate, heavier validation this script does not attempt.

Writes a timestamped JSON snapshot to testsys/parity/evidence_output/
(gitignored), like the TPV29/TPV34 scripts.
"""
import argparse
import datetime
import glob
import json
import math
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TESTSYS = os.path.dirname(HERE)
REPO_ROOT = os.path.dirname(TESTSYS)
OUT_DIR = os.path.join(HERE, 'evidence_output')
SUBMISSION_DIR_DEFAULT = os.path.join(REPO_ROOT, 'scec_archive', 'tpv35',
                                      'eqdyna3d-v4.1.3-100m-2017')
STATIONS_FILE_DEFAULT = os.path.join(REPO_ROOT, 'case_input', 'test.tpv35',
                                     'tpv35_station_locations.txt')
SLIP_RATE_THRESHOLD = 1.0e-3   # m/s, SCEC rupture-front definition
UNRUPTURED_CPLOT = 1000.0      # s, the 2017 cplot sentinel
UNRUPTURED_FRT = 99999.0       # s, EQdyna frt column 4 sentinel
FAULT_AREA_KM2 = 40.0 * 15.5
ON_FAULT_X_KM = (-25.0, -20.0, -15.0, -10.0, -5.0, 0.0, 5.0)
ON_FAULT_DEPTH_KM_SPEC = 8.1
OFF_FAULT_RATIO_FLOOR = 0.01


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
    """Station-name field in hectometres, as the Fortran writer formats it
    (library_output.f90's i4.3 edit: >=3 digits, sign carried separately)."""
    v = int(round(km * 10.0))
    return '%03d' % v if v >= 0 else '-%03d' % (-v)


def _nint(v):
    """Fortran NINT: round half away from zero (Python round() is banker's
    rounding and disagrees at exact .5 ties)."""
    return int(math.floor(v + 0.5)) if v >= 0 else int(math.ceil(v - 0.5))


def onFaultStationPairs(gate_dz_m):
    """[(run_name, archive_name, x_km, depth_km_gate), ...] for the 7
    on-fault stations. The archive is itself a 100 m EQdyna run, so its
    depth label is always dp081 (8.1 km); the run's own label depends on
    ITS dz (dp081 at 100 m, dp080 at the 500 m gate config)."""
    dg = round(ON_FAULT_DEPTH_KM_SPEC * 1.0e3 / gate_dz_m) * gate_dz_m / 1.0e3
    out = []
    for x in ON_FAULT_X_KM:
        out.append(('faultst%sdp%s.txt' % (_hm(x), _hm(dg)),
                    'faultst%sdp%s' % (_hm(x), _hm(ON_FAULT_DEPTH_KM_SPEC)),
                    x, dg))
    return out


def readOffFaultStations(path=STATIONS_FILE_DEFAULT):
    """[(name, x_eq_km, y_eq_km), ...] in file order. x_eq = x_tpv (along
    strike), y_eq = z_tpv (fault-normal) -- tpv35Tools.py's frame mapping,
    confirmed against TPV35_Description_v05.pdf Part 2 text (near side
    z<0, far side z>0, matching tpv35Tools.py's y_eq<0/y_eq>0 exactly)."""
    out = []
    with open(path) as f:
        for line in f:
            s = line.split()
            if not s or s[0] == 'station_name':
                continue
            z_tpv, x_tpv = float(s[2]), float(s[3])
            out.append((s[0], x_tpv / 1.0e3, z_tpv / 1.0e3))
    return out


def offFaultBodyFilename(x_eq_km, y_eq_km):
    """The run's own off-fault file name for a station at this (x_eq, y_eq),
    z_eq = 0 (surface) -- library_output.f90:217-221: 'body' <- y4nds (hm),
    'st' <- x4nds (hm), 'dp' <- |z4nds| (hm), from the REQUESTED coordinate,
    not the snapped node (row 94 ruling), so this is resolution-independent."""
    st_hm = _nint(x_eq_km * 10.0)
    body_hm = _nint(y_eq_km * 10.0)
    def fmt3(v):
        return '%03d' % v if v >= 0 else '-%03d' % (-v)
    return 'body%sst%sdp000.txt' % (fmt3(body_hm), fmt3(st_hm))


def first_sample_stress(path):
    a = _numeric_rows(path)
    return float(a[0, 3]), float(a[0, 7]), float(a[0, 0])


def arrival_time(path, term):
    a = _numeric_rows(path)
    a = a[a[:, 0] <= term + 1e-9]
    hit = np.nonzero(a[:, 2] > SLIP_RATE_THRESHOLD)[0]
    return float(a[hit[0], 0]) if hit.size else None


def slip_at(path, t):
    a = _numeric_rows(path)
    i = int(np.argmin(np.abs(a[:, 0] - t)))
    return float(a[i, 1]), float(a[i, 0])


def peak_vel(path, term, hcol=2, vcol=4):
    """(peak |h-vel|, t of peak, peak |v-vel|, t of peak) within term, from
    an off-fault body/named-station file (columns t h-disp h-vel v-disp
    v-vel n-disp n-vel)."""
    a = _numeric_rows(path)
    a = a[a[:, 0] <= term + 1e-9]
    ih = int(np.argmax(np.abs(a[:, hcol])))
    iv = int(np.argmax(np.abs(a[:, vcol])))
    return (float(abs(a[ih, hcol])), float(a[ih, 0]),
            float(abs(a[iv, vcol])), float(a[iv, 0]))


def ruptured_fraction_cplot(path, term):
    a = _numeric_rows(path)
    t = a[:, 2]
    return float(np.mean((t < UNRUPTURED_CPLOT) & (t <= term))), int(a.shape[0])


def ruptured_fraction_frt(run_dir, term):
    files = sorted(glob.glob(os.path.join(run_dir, 'frt.txt*')))
    if not files:
        raise SystemExit('%s has no frt.txt* files' % run_dir)
    a = np.vstack([np.loadtxt(f) for f in files])
    _, first = np.unique(np.round(a[:, :3], 3), axis=0, return_index=True)
    t = a[first, 3]
    return float(np.mean((t < UNRUPTURED_FRT) & (t <= term))), int(first.size)


def run_term_and_dz(run_dir):
    """Read the run's term (s) and dz (m) off an on-fault station header --
    globs faultst000dp*.txt since the depth label depends on the run's own
    dz (unknown until read)."""
    cands = sorted(glob.glob(os.path.join(run_dir, 'faultst000dp*.txt')))
    if not cands:
        raise SystemExit('%s has no faultst000dp*.txt to read term/dz from' % run_dir)
    path = cands[0]
    dz = term = dt = None
    with open(path) as f:
        for line in f:
            key = line.lstrip('# ').split('=')[0].strip()
            if key == 'element_size':
                dz = float(line.split('=')[1])
            elif key == 'time_step':
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


def compare(sub, run, stations_file):
    term, dz = run_term_and_dz(run)
    r = dict(run_term_s=term, run_dz_m=dz, on_fault=[], off_fault=[])

    for rname, aname, x, dg in onFaultStationPairs(dz):
        rp, ap = os.path.join(run, rname), os.path.join(sub, aname)
        if not (os.path.isfile(rp) and os.path.isfile(ap)):
            r['on_fault'].append(dict(run=rname, archive=aname, missing=True))
            continue
        gs, gn, gt0 = first_sample_stress(rp)
        as_, an, at0 = first_sample_stress(ap)
        r['on_fault'].append(dict(
            run=rname, archive=aname, x_km=x, depth_km=dg,
            shear0_run_MPa=gs, shear0_2017_MPa=as_, shear0_ratio=gs / as_,
            normal0_run_MPa=gn, normal0_2017_MPa=an, normal0_ratio=gn / an,
            first_sample_s=(gt0, at0),
            arrival_run_s=arrival_time(rp, term), arrival_2017_s=arrival_time(ap, term)))
    ok = [s for s in r['on_fault'] if not s.get('missing')]
    sr = np.array([s['shear0_ratio'] for s in ok])
    nr = np.array([s['normal0_ratio'] for s in ok])
    r['initial_stress'] = dict(
        n=len(ok), shear_ratio_median=float(np.median(sr)), shear_ratio_min=float(sr.min()),
        shear_ratio_max=float(sr.max()), normal_ratio_median=float(np.median(nr)),
        normal_ratio_min=float(nr.min()), normal_ratio_max=float(nr.max()),
        sign_agreement=bool(np.all(sr > 0) and np.all(nr > 0)))
    both = [s for s in ok if s['arrival_run_s'] is not None and s['arrival_2017_s'] is not None]
    dts = np.array([s['arrival_run_s'] - s['arrival_2017_s'] for s in both])
    r['arrival'] = dict(
        n_both=len(both), n_run_only=sum(1 for s in ok if s['arrival_run_s'] is not None
                                         and s['arrival_2017_s'] is None),
        n_2017_only=sum(1 for s in ok if s['arrival_2017_s'] is not None
                        and s['arrival_run_s'] is None),
        median_dt_s=float(np.median(dts)) if dts.size else None,
        max_abs_dt_s=float(np.abs(dts).max()) if dts.size else None,
        median_abs_dt_s=float(np.median(np.abs(dts))) if dts.size else None)

    fc, nc = ruptured_fraction_cplot(os.path.join(sub, 'cplot'), term)
    ff, nf = ruptured_fraction_frt(run, term)
    r['ruptured_fraction_at_term'] = dict(archive_2017=fc, archive_nodes=nc, run=ff, run_nodes=nf)

    hyp_run = os.path.join(run, [n for n, _, _, _ in onFaultStationPairs(dz) if 'st000' in n][0])
    r['hypocentre_slip_at_term_m'] = dict(
        run=slip_at(hyp_run, term), archive_2017=slip_at(os.path.join(sub, 'faultst000dp081'), term))

    # Two requested stations can round to the SAME body file name (4074_viney
    # and 4134_vyc are 54 m apart: both are body-039st-183): the solver keeps
    # the first and drops the second with a warning, so the file belongs to
    # the FIRST station only. Never credit it to the second.
    owner = {}
    for name, x_eq, y_eq in readOffFaultStations(stations_file):
        ap = os.path.join(sub, name)
        rfile = offFaultBodyFilename(x_eq, y_eq)
        rp = os.path.join(run, rfile)
        collides = owner.setdefault(rfile, name) != name
        run_present = os.path.isfile(rp) and not collides
        if not (os.path.isfile(ap) and run_present):
            r['off_fault'].append(dict(name=name, run_file=rfile,
                                       archive_present=os.path.isfile(ap),
                                       run_present=run_present, missing=True,
                                       run_file_belongs_to=(owner[rfile] if collides else None)))
            continue
        hr, thr, vr, tvr = peak_vel(rp, term)
        ha, tha, va, tva = peak_vel(ap, term)
        r['off_fault'].append(dict(
            name=name, run_file=os.path.basename(rp), x_eq_km=x_eq, y_eq_km=y_eq,
            peak_hvel_run=hr, peak_hvel_2017=ha, peak_hvel_ratio=(hr / ha if ha else None),
            peak_hvel_dt_s=thr - tha,
            peak_vvel_run=vr, peak_vvel_2017=va, peak_vvel_ratio=(vr / va if va else None),
            peak_vvel_dt_s=tvr - tva))
    ofok = [s for s in r['off_fault'] if not s.get('missing')]
    # A ratio is only meaningful where the 2017 station has actually been
    # reached by the wavefield inside the term: below 1% of the largest 2017
    # peak over all stations, the denominator is numerical noise (a short-
    # term run gave ratios of 1e13 there). Those stations stay in the JSON
    # with their raw peaks; they are just left out of the ratio statistics.
    hmax = max([s['peak_hvel_2017'] for s in ofok] or [0.0])
    vmax = max([s['peak_vvel_2017'] for s in ofok] or [0.0])
    for s in ofok:
        s['hvel_in_stats'] = s['peak_hvel_2017'] >= OFF_FAULT_RATIO_FLOOR * hmax > 0
        s['vvel_in_stats'] = s['peak_vvel_2017'] >= OFF_FAULT_RATIO_FLOOR * vmax > 0
    hvr = np.array([s['peak_hvel_ratio'] for s in ofok if s['hvel_in_stats']])
    vvr = np.array([s['peak_vvel_ratio'] for s in ofok if s['vvel_in_stats']])
    hdt = np.array([s['peak_hvel_dt_s'] for s in ofok if s['hvel_in_stats']])
    vdt = np.array([s['peak_vvel_dt_s'] for s in ofok if s['vvel_in_stats']])
    r['off_fault_summary'] = dict(
        n=len(ofok), n_missing=len(r['off_fault']) - len(ofok),
        ratio_floor_fraction_of_max=OFF_FAULT_RATIO_FLOOR,
        n_hvel_in_stats=int(hvr.size), n_vvel_in_stats=int(vvr.size),
        hvel_peak_time_median_abs_dt_s=float(np.median(np.abs(hdt))) if hdt.size else None,
        vvel_peak_time_median_abs_dt_s=float(np.median(np.abs(vdt))) if vdt.size else None,
        hvel_ratio_median=float(np.median(hvr)) if hvr.size else None,
        hvel_ratio_min=float(hvr.min()) if hvr.size else None,
        hvel_ratio_max=float(hvr.max()) if hvr.size else None,
        vvel_ratio_median=float(np.median(vvr)) if vvr.size else None,
        vvel_ratio_min=float(vvr.min()) if vvr.size else None,
        vvel_ratio_max=float(vvr.max()) if vvr.size else None)
    return r


def print_report(r, sub, run):
    print('TPV35: run %s (dz %g m, term %g s) vs 2017 SCEC submission %s (100 m, 18 s)'
          % (run, r['run_dz_m'], r['run_term_s'], sub))
    i = r['initial_stress']
    print('\n1. initial stress, run/2017 ratio over %d on-fault stations:' % i['n'])
    print('   h-shear  median %.4f  range %.4f .. %.4f' % (i['shear_ratio_median'],
                                                           i['shear_ratio_min'], i['shear_ratio_max']))
    print('   n-stress median %.4f  range %.4f .. %.4f   signs agree: %s'
          % (i['normal_ratio_median'], i['normal_ratio_min'], i['normal_ratio_max'],
             i['sign_agreement']))
    print('   %-20s %-18s %9s %9s %9s %9s' % ('run', '2017', 'tau_r', 'tau_17', 'sn_r', 'sn_17'))
    for s in r['on_fault']:
        if s.get('missing'):
            print('   %-20s %-18s MISSING' % (s['run'], s['archive']))
            continue
        print('   %-20s %-18s %9.3f %9.3f %9.3f %9.3f' % (
            s['run'], s['archive'], s['shear0_run_MPa'], s['shear0_2017_MPa'],
            s['normal0_run_MPa'], s['normal0_2017_MPa']))
    a = r['arrival']
    print('\n2. rupture arrival (h-slip-rate > 1 mm/s) within %g s: %d stations in both, '
          '%d run-only, %d 2017-only' % (r['run_term_s'], a['n_both'], a['n_run_only'],
                                         a['n_2017_only']))
    if a['n_both']:
        print('   run - 2017: median %+.3f s, median |dt| %.3f s, max |dt| %.3f s'
              % (a['median_dt_s'], a['median_abs_dt_s'], a['max_abs_dt_s']))
        for s in r['on_fault']:
            if not s.get('missing') and s['arrival_run_s'] is not None and s['arrival_2017_s'] is not None:
                print('   %-20s run %6.3f  2017 %6.3f  dt %+.3f' % (
                    s['run'], s['arrival_run_s'], s['arrival_2017_s'],
                    s['arrival_run_s'] - s['arrival_2017_s']))
    f = r['ruptured_fraction_at_term']
    print('\n3. fraction of fault ruptured by %g s: 2017 cplot %.3f (%d nodes), run frt %.3f (%d nodes)'
          % (r['run_term_s'], f['archive_2017'], f['archive_nodes'], f['run'], f['run_nodes']))
    h = r['hypocentre_slip_at_term_m']
    print('4. hypocentre h-slip at %g s: run %.3f m (t=%.3f), 2017 %.3f m (t=%.3f)'
          % (r['run_term_s'], h['run'][0], h['run'][1], h['archive_2017'][0], h['archive_2017'][1]))
    o = r['off_fault_summary']
    print('\n5. off-fault peak ground velocity, run/2017 ratio over %d of %d named stations '
          '(%d missing one side or the other):' % (o['n'], len(r['off_fault']), o['n_missing']))
    print('   (ratio stats only where the 2017 peak >= %g x the largest 2017 peak)'
          % o['ratio_floor_fraction_of_max'])
    if o['n_hvel_in_stats']:
        print('   peak |h-vel| over %2d: median %.4f  range %.4f .. %.4f  median |t_peak dt| %.3f s'
              % (o['n_hvel_in_stats'], o['hvel_ratio_median'], o['hvel_ratio_min'],
                 o['hvel_ratio_max'], o['hvel_peak_time_median_abs_dt_s']))
    if o['n_vvel_in_stats']:
        print('   peak |v-vel| over %2d: median %.4f  range %.4f .. %.4f  median |t_peak dt| %.3f s'
              % (o['n_vvel_in_stats'], o['vvel_ratio_median'], o['vvel_ratio_min'],
                 o['vvel_ratio_max'], o['vvel_peak_time_median_abs_dt_s']))
    print('   %-12s %-24s %9s %9s %7s %9s %9s %7s' % ('station', 'run file', 'hv_run', 'hv_17',
                                                   'ratio', 'vv_run', 'vv_17', 'ratio'))
    for s in r['off_fault']:
        if not s.get('missing'):
            print('   %-12s %-24s %9.4f %9.4f %7.3f %9.4f %9.4f %7.3f'
                  % (s['name'], s['run_file'], s['peak_hvel_run'], s['peak_hvel_2017'],
                     s['peak_hvel_ratio'] or float('nan'), s['peak_vvel_run'],
                     s['peak_vvel_2017'], s['peak_vvel_ratio'] or float('nan')))
    for s in r['off_fault']:
        if s.get('missing'):
            print('   %-12s run_file=%-22s archive_present=%s run_present=%s MISSING%s'
                  % (s['name'], s['run_file'], s['archive_present'], s['run_present'],
                     ('  (name collides with %s, whose file it is)' % s['run_file_belongs_to'])
                     if s.get('run_file_belongs_to') else ''))
    print('\nThis is a report (NOT a comparison to real Parkfield recordings -- see module '
          'docstring), not a gate.')


OFF_FAULT_TS_N = 6   # off-fault ts overlay: the N stations with the largest 2017 peak |h-vel|


def _frt_points(run_dir):
    files = sorted(glob.glob(os.path.join(run_dir, 'frt.txt*')))
    a = np.vstack([np.loadtxt(f, usecols=(0, 1, 2, 3)) for f in files])
    _, first = np.unique(np.round(a[:, :3], 3), axis=0, return_index=True)
    return a[first]


def make_overlays(r, sub, run, plot_dir):
    """Owner 2026-10-06 (board row 19(c)): cplot and ts overlays, run vs 2017.
    Writes three PNGs into plot_dir and returns their absolute paths."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.tri
    from matplotlib.lines import Line2D
    os.makedirs(plot_dir, exist_ok=True)
    term, dz = r['run_term_s'], r['run_dz_m']
    out = []

    # (1) cplot: rupture-time contours on the fault, same axes, same levels.
    c = _numeric_rows(os.path.join(sub, 'cplot'))
    xs, ks = np.unique(c[:, 0]), np.unique(c[:, 1])
    T = np.full((ks.size, xs.size), np.nan)
    T[np.searchsorted(ks, c[:, 1]), np.searchsorted(xs, c[:, 0])] = c[:, 2]
    T[T >= UNRUPTURED_CPLOT] = np.nan
    f = _frt_points(run)
    m = f[:, 3] < UNRUPTURED_FRT
    levels = np.arange(0.5, term + 1e-9, 0.5)
    fig, ax = plt.subplots(figsize=(12, 4.2))
    ax.contour(xs / 1e3, ks / 1e3, T, levels=levels, colors='k', linewidths=0.9)
    # Delaunay triangulation of only the RUPTURED points fills the convex hull
    # of that (non-convex) patch, so a handful of long thin triangles bridge
    # across unruptured bays -- tricontour then draws spurious lines through
    # them (the diagonal streaks the owner flagged, 2026-10-06). Mask any
    # triangle whose longest edge exceeds a few node spacings; 2*dz is wide
    # enough for this patch's own node spacing but rejects a bridge across
    # unruptured ground.
    tri = matplotlib.tri.Triangulation(f[m, 0] / 1e3, -f[m, 2] / 1e3)
    pts = np.column_stack([tri.x, tri.y])
    edges = pts[tri.triangles] - pts[tri.triangles[:, [1, 2, 0]]]
    max_edge_km = np.max(np.hypot(edges[:, :, 0], edges[:, :, 1]), axis=1)
    tri.set_mask(max_edge_km > 2.0 * dz / 1e3)
    ax.tricontour(tri, f[m, 3], levels=levels, colors='r', linewidths=0.9, linestyles='--')
    for x in ON_FAULT_X_KM:
        ax.plot(x, ON_FAULT_DEPTH_KM_SPEC, 'b^', ms=6)
    ax.invert_yaxis()
    ax.set_aspect('equal')
    ax.set_xlabel('along strike (km)')
    ax.set_ylabel('down dip (km)')
    ax.set_title('TPV35 rupture-time contours every 0.5 s to %.1f s: 2017 100 m (black) vs run %g m '
                 '(red dashed); ruptured fraction 2017 %.3f, run %.3f'
                 % (term, dz, r['ruptured_fraction_at_term']['archive_2017'],
                    r['ruptured_fraction_at_term']['run']), fontsize=9)
    ax.legend([Line2D([], [], color='k'), Line2D([], [], color='r', ls='--'),
               Line2D([], [], color='b', marker='^', ls='')],
              ['2017 submission', 'this run', 'on-fault stations'], fontsize=8, loc='lower left')
    pth = os.path.abspath(os.path.join(plot_dir, 'tpv35_cplot_overlay.png'))
    fig.tight_layout(); fig.savefig(pth, dpi=150); plt.close(fig); out.append(pth)

    # (2) on-fault ts: h-slip-rate (col 2) and h-shear-stress (col 3), all 7 stations.
    pairs = onFaultStationPairs(dz)
    fig, axs = plt.subplots(len(pairs), 2, figsize=(11, 2.0 * len(pairs)), sharex=True)
    for i, (rn, an, xkm, _) in enumerate(pairs):
        a = _numeric_rows(os.path.join(sub, an)); b = _numeric_rows(os.path.join(run, rn))
        for j, (col, lab) in enumerate(((2, 'h-slip-rate (m/s)'), (3, 'h-shear-stress (MPa)'))):
            axs[i, j].plot(a[:, 0], a[:, col], 'k', lw=1, label='2017')
            axs[i, j].plot(b[:, 0], b[:, col], 'r--', lw=1, label='run')
            axs[i, j].set_ylabel(lab, fontsize=8)
            axs[i, j].set_title('%s  (x=%g km)' % (an, xkm), fontsize=8)
    axs[0, 0].legend(fontsize=8)
    for a_ in axs[-1]:
        a_.set_xlabel('t (s)'); a_.set_xlim(0, term)
    fig.tight_layout()
    pth = os.path.abspath(os.path.join(plot_dir, 'tpv35_ts_onfault_overlay.png'))
    fig.savefig(pth, dpi=120); plt.close(fig); out.append(pth)

    # (3) off-fault ts: h-vel (col 2) and v-vel (col 4) at the OFF_FAULT_TS_N strongest stations.
    ok = [s_ for s_ in r['off_fault'] if not s_.get('missing')]
    ok = sorted(ok, key=lambda s_: -s_['peak_hvel_2017'])[:OFF_FAULT_TS_N]
    fig, axs = plt.subplots(len(ok), 2, figsize=(11, 2.0 * len(ok)), sharex=True)
    for i, s_ in enumerate(ok):
        a = _numeric_rows(os.path.join(sub, s_['name'])); b = _numeric_rows(os.path.join(run, s_['run_file']))
        for j, (col, lab) in enumerate(((2, 'h-vel (m/s)'), (4, 'v-vel (m/s)'))):
            axs[i, j].plot(a[:, 0], a[:, col], 'k', lw=1, label='2017')
            axs[i, j].plot(b[:, 0], b[:, col], 'r--', lw=1, label='run')
            axs[i, j].set_ylabel(lab, fontsize=8)
            axs[i, j].set_title('%s  (x=%.1f km, normal %.1f km)' % (s_['name'], s_['x_eq_km'],
                                s_['y_eq_km']), fontsize=8)
    axs[0, 0].legend(fontsize=8)
    for a_ in axs[-1]:
        a_.set_xlabel('t (s)'); a_.set_xlim(0, term)
    fig.tight_layout()
    pth = os.path.abspath(os.path.join(plot_dir, 'tpv35_ts_offfault_overlay.png'))
    fig.savefig(pth, dpi=120); plt.close(fig); out.append(pth)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--submission-dir', default=SUBMISSION_DIR_DEFAULT)
    ap.add_argument('--stations-file', default=STATIONS_FILE_DEFAULT)
    ap.add_argument('--run-dir', required=True,
                    help='a COMPLETED test.tpv35 run directory at the SPEC resolution/term '
                         '(dx=100, term=18 -- edit a scratch case copy of '
                         'user_defined_params.py, never the committed gate config): '
                         'create.newcase <dir> test.tpv35; cd <dir>; edit par.dx/par.term; '
                         'python3 case.setup; mpirun -np N eqdyna')
    ap.add_argument('--no-json', action='store_true')
    ap.add_argument('--plot-dir', default=None,
                    help='write the cplot and ts overlay PNGs here (owner 2026-10-06, row 19(c))')
    args = ap.parse_args()
    if not os.path.isdir(args.submission_dir):
        raise SystemExit('%s missing -- fetch per scec_archive/tpv35/*/PROVENANCE.md'
                         % args.submission_dir)
    r = compare(args.submission_dir, args.run_dir, args.stations_file)
    print_report(r, args.submission_dir, args.run_dir)
    if args.plot_dir:
        r['overlay_pngs'] = make_overlays(r, args.submission_dir, args.run_dir, args.plot_dir)
        for pth in r['overlay_pngs']:
            print('overlay: %s' % pth)
    if not args.no_json:
        os.makedirs(OUT_DIR, exist_ok=True)
        stamp = datetime.datetime.utcnow().strftime('%Y%m%dT%H%M%SZ')
        path = os.path.join(OUT_DIR, 'tpv35_scec_comparison_%s.json' % stamp)
        with open(path, 'w') as f:
            json.dump(dict(provenance=provenance(), submission_dir=args.submission_dir,
                           run_dir=args.run_dir, result=r), f, indent=1)
        print('snapshot: %s' % path)
    return 0


if __name__ == '__main__':
    sys.exit(main())

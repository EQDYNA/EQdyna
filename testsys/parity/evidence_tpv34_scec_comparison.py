#! /usr/bin/env python3
"""
Compares a test.tpv34 GATE-configuration run (500 m, 5 s, 4 ranks) with
EQdyna's own 2016 SCEC TPV34 submission (50 m, 20 s), FROM THE ARCHIVED
DATA -- the rule 17 step 6 independent check for test.tpv34, committed as a
script so its numbers carry provenance instead of living in prose
(PROJECT_RULES rule 4 / rule 6; the TPV29 lesson, pathway item 28).

REPORT-ONLY, same as evidence_tpv29_scec_comparison.py: never a gate, never
wired into testsys/run.py, invoked by hand. Prints the numbers, says what
each one tests and what it cannot test, writes a JSON snapshot, and always
exits 0.

WHAT IS COMPARED
----------------
The two runs differ by a factor 10 in dx and 4 in duration, and the 2016
run used CVM-H 15.1.0 where the gate's shipped grid is 15.1.1 (see
case_input/test.tpv34/README.md), so only what is comparable IS compared:

  1. INITIAL STRESS at the 35 on-fault stations: first-sample h-shear-stress
     and n-stress, gate run vs 2016 (the 2016 files start at t = 0.004 s,
     before any slip at every station but the hypocentre's immediate
     neighbours; the gate's first sample is at its own dt). This tests the
     whole Part 2 + Part 3 chain -- CVM-H extraction, UTM mapping, clamp,
     mu/mu0 scaling, nucleation increment -- against an independent 2016
     extraction at 10x finer sampling. Expected: agreement to a few percent,
     with the gate's 8-cell 500 m average differing from a 50 m average most
     where the medium varies fastest (shallow sediments). A sign flip or a
     tens-of-percent miss would mean a wrong frame, wrong UTM, or wrong
     scaling.
  2. RUPTURE ARRIVAL at the on-fault stations that rupture within the gate
     term in BOTH runs: first time h-slip-rate exceeds SCEC's 1 mm/s
     threshold, t_gate - t_2016, listed and summarised. Mesh and medium
     differ (500 m vs 50 m; 15.1.1 vs 15.1.0), so this is a
     direction-and-magnitude report, not an accuracy claim. Measured
     2026-10-06: the 500 m front arrives LATER everywhere, median +0.30 s,
     max +0.54 s at 16 stations reached within 5 s in both runs.
  3. RUPTURED AREA AT THE GATE TERM: fraction of the 30 x 15 km fault with
     rupture time <= term from the 2016 cplot (50 m nodes) vs the gate run's
     frt files (500 m nodes), both node-counted on their own grids.
  4. HYPOCENTRE SLIP AT THE GATE TERM: h-slip at t = term at faultst000dp075.

The 2016 archive's 2.4 km stations are named dp024; the gate's are dp025
(500 m fault planes). They are paired by the gate's rounding rule and the
pairing is printed; the 500 m station is 100 m deeper than the 2016 one.

Writes a timestamped JSON snapshot to testsys/parity/evidence_output/
(gitignored), like the TPV29 script.
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


def station_pairs(gate_dz_m):
    """[(gate_name, archive_name, x_km, depth_km_spec, depth_km_gate), ...]"""
    out = []
    for x in ON_FAULT_X_KM:
        for d in ON_FAULT_DEPTH_KM_SPEC:
            dg = round(d * 1.0e3 / gate_dz_m) * gate_dz_m / 1.0e3
            out.append(('faultst%sdp%s.txt' % (_hm(x), _hm(dg)),
                        'faultst%sdp%s' % (_hm(x), _hm(d)), x, d, dg))
    return out


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


def gate_term_and_dz(run_dir):
    """Read the run's term (s) and dz (m) off a station file header."""
    path = os.path.join(run_dir, 'faultst000dp075.txt')
    dz = term = None
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


def compare(sub, run):
    term, dz = gate_term_and_dz(run)
    r = dict(gate_term_s=term, gate_dz_m=dz, stations=[])
    for gname, aname, x, d, dg in station_pairs(dz):
        gp, apath = os.path.join(run, gname), os.path.join(sub, aname)
        if not (os.path.isfile(gp) and os.path.isfile(apath)):
            r['stations'].append(dict(gate=gname, archive=aname, missing=True))
            continue
        gs, gn, gt0 = first_sample_stress(gp)
        as_, an, at0 = first_sample_stress(apath)
        r['stations'].append(dict(
            gate=gname, archive=aname, x_km=x, depth_km_spec=d, depth_km_gate=dg,
            shear0_gate_MPa=gs, shear0_2016_MPa=as_, shear0_ratio=gs / as_,
            normal0_gate_MPa=gn, normal0_2016_MPa=an, normal0_ratio=gn / an,
            first_sample_s=(gt0, at0),
            arrival_gate_s=arrival_time(gp, term), arrival_2016_s=arrival_time(apath, term)))
    ok = [s for s in r['stations'] if not s.get('missing')]
    sr = np.array([s['shear0_ratio'] for s in ok])
    nr = np.array([s['normal0_ratio'] for s in ok])
    r['initial_stress'] = dict(
        n=len(ok), shear_ratio_median=float(np.median(sr)), shear_ratio_min=float(sr.min()),
        shear_ratio_max=float(sr.max()), normal_ratio_median=float(np.median(nr)),
        normal_ratio_min=float(nr.min()), normal_ratio_max=float(nr.max()),
        sign_agreement=bool(np.all(sr > 0) and np.all(nr > 0)))
    both = [s for s in ok if s['arrival_gate_s'] is not None and s['arrival_2016_s'] is not None]
    dts = np.array([s['arrival_gate_s'] - s['arrival_2016_s'] for s in both])
    r['arrival'] = dict(
        n_both=len(both), n_gate_only=sum(1 for s in ok if s['arrival_gate_s'] is not None
                                          and s['arrival_2016_s'] is None),
        n_2016_only=sum(1 for s in ok if s['arrival_2016_s'] is not None
                        and s['arrival_gate_s'] is None),
        median_dt_s=float(np.median(dts)) if dts.size else None,
        max_abs_dt_s=float(np.abs(dts).max()) if dts.size else None,
        median_abs_dt_s=float(np.median(np.abs(dts))) if dts.size else None)
    fc, nc = ruptured_fraction_cplot(os.path.join(sub, 'cplot'), term)
    ff, nf = ruptured_fraction_frt(run, term)
    r['ruptured_fraction_at_term'] = dict(archive_2016=fc, archive_nodes=nc, gate=ff, gate_nodes=nf)
    r['hypocentre_slip_at_term_m'] = dict(
        gate=slip_at(os.path.join(run, 'faultst000dp075.txt'), term),
        archive_2016=slip_at(os.path.join(sub, 'faultst000dp075'), term))
    return r


def print_report(r, sub, run):
    print('TPV34: gate run %s (dx %g m, term %g s) vs 2016 SCEC submission %s (50 m, 20 s)'
          % (run, r['gate_dz_m'], r['gate_term_s'], sub))
    i = r['initial_stress']
    print('\n1. initial stress, gate/2016 ratio over %d on-fault stations:' % i['n'])
    print('   h-shear  median %.4f  range %.4f .. %.4f' % (i['shear_ratio_median'],
                                                           i['shear_ratio_min'], i['shear_ratio_max']))
    print('   n-stress median %.4f  range %.4f .. %.4f   signs agree: %s'
          % (i['normal_ratio_median'], i['normal_ratio_min'], i['normal_ratio_max'],
             i['sign_agreement']))
    print('   %-22s %-18s %9s %9s %9s %9s' % ('gate', '2016', 'tau_g', 'tau_16', 'sn_g', 'sn_16'))
    for s in r['stations']:
        if s.get('missing'):
            print('   %-22s %-18s MISSING' % (s['gate'], s['archive']))
            continue
        print('   %-22s %-18s %9.3f %9.3f %9.3f %9.3f' % (
            s['gate'], s['archive'], s['shear0_gate_MPa'], s['shear0_2016_MPa'],
            s['normal0_gate_MPa'], s['normal0_2016_MPa']))
    a = r['arrival']
    print('\n2. rupture arrival (h-slip-rate > 1 mm/s) within %g s: %d stations in both, '
          '%d gate-only, %d 2016-only' % (r['gate_term_s'], a['n_both'], a['n_gate_only'],
                                          a['n_2016_only']))
    if a['n_both']:
        print('   gate - 2016: median %+.3f s, median |dt| %.3f s, max |dt| %.3f s'
              % (a['median_dt_s'], a['median_abs_dt_s'], a['max_abs_dt_s']))
        for s in r['stations']:
            if not s.get('missing') and s['arrival_gate_s'] is not None and s['arrival_2016_s'] is not None:
                print('   %-22s gate %6.3f  2016 %6.3f  dt %+.3f' % (
                    s['gate'], s['arrival_gate_s'], s['arrival_2016_s'],
                    s['arrival_gate_s'] - s['arrival_2016_s']))
    f = r['ruptured_fraction_at_term']
    print('\n3. fraction of fault ruptured by %g s: 2016 cplot %.3f (%d nodes), gate frt %.3f (%d nodes)'
          % (r['gate_term_s'], f['archive_2016'], f['archive_nodes'], f['gate'], f['gate_nodes']))
    h = r['hypocentre_slip_at_term_m']
    print('4. hypocentre h-slip at %g s: gate %.3f m (t=%.3f), 2016 %.3f m (t=%.3f)'
          % (r['gate_term_s'], h['gate'][0], h['gate'][1], h['archive_2016'][0], h['archive_2016'][1]))
    print('\nThis is a report (500 m / 5 s vs 50 m / 20 s, CVM-H 15.1.1 vs 15.1.0), not a gate.')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--submission-dir', default=SUBMISSION_DIR_DEFAULT)
    ap.add_argument('--run-dir', required=True,
                    help='a COMPLETED test.tpv34 gate-config run directory (faultst*.txt, '
                         'frt.txt*): create.newcase <dir> test.tpv34; cd <dir>; python3 '
                         'case.setup; mpirun -np 4 eqdyna')
    ap.add_argument('--no-json', action='store_true')
    args = ap.parse_args()
    if not os.path.isdir(args.submission_dir):
        raise SystemExit('%s missing -- fetch per scec_archive/tpv34/*/PROVENANCE.md'
                         % args.submission_dir)
    r = compare(args.submission_dir, args.run_dir)
    print_report(r, args.submission_dir, args.run_dir)
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

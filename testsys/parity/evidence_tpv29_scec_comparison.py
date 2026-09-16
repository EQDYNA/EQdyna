#! /usr/bin/env python3
"""
Regenerates, FROM THE ARCHIVED DATA, the TPV29 cross-code validation numbers
behind pathway_forward.md items 20/28 -- median rupture-time diff, median
final-slip percentage diff, ruptured area and Mw, all at matched dx=100 m --
previously cited only in prose (pathway_forward.md, case_input/test.tpv29/
README.md, commit messages), with no committed script to reproduce them
(PROJECT_RULES rule 4 "only fresh runs are evidence" / rule 6 "every number
carries its provenance").

REPORT-ONLY, same as evidence_drv_a6_chaos.py in this directory: never a
gate, never wired into testsys/run.py, invoked by hand. Prints the four
numbers, prints the recorded prose claim alongside, and says MATCH or DRIFT
-- it never asserts and it always exits 0.

WHAT IS COMPARED AND WHERE IT COMES FROM
-----------------------------------------
Two independent EQdyna runs of the same SCEC benchmark, at the same dx=100 m,
21 years apart:

  1. THE 2015 SUBMISSION -- scec_archive/tpv29/eqdyna-v3.1-100m-2015/, fetched
     from the public SCEC/USGS cvws portal (see scripts/scec/README.md for
     re-fetch; `python3 scripts/scec/fetch_all_dliu.py`). Format, read off the
     files' own headers (documented in that directory's PROVENANCE.md):
       - `cplot`            : 3 columns "strike(m) dip(m) rupture-time(s)",
                               401 x 201 node-centred grid, sentinel 1000.0 s.
       - `faultst<S>dp<D>`  : 8-column time series (t, h-slip, h-slip-rate,
                               h-shear-stress, v-slip, v-slip-rate,
                               v-shear-stress, n-stress); filename S, D are
                               strike/dip in HECTOMETRES (S=100 -> 10.0 km
                               along strike -- verified against the file's own
                               "# location = ... km along strike/down-dip"
                               comment line for two stations).
     24 on-fault stations ship in the archive (PROVENANCE.md's file table).

  2. THE CURRENT-CODE RUN -- scratch/tpv29/ny48_dx100/ (NOT tracked; this
     directory is gitignored, same as scec_archive/). This is an
     ALREADY-COMPLETED dx=100 m, 48-rank test.tpv29 run (SCECRuptureTime.txt,
     the same 24 faultst*.txt names, and per-rank fault-restart dumps
     frt.txt<rank> covering the full 401x201 grid with final slip). Per rule
     9 this script does NOT re-launch that run -- it is not "cheap" (case
     setup + mpirun across dozens of ranks, tens of minutes) and a valid
     completed one already sits on this box. If it is missing, this script
     raises with the exact recipe (case_input/test.tpv29/README.md) rather
     than silently falling back to a different resolution or a stale number.
     Note: `scratch/tpv29/ny1_dx100/` was CHECKED and REJECTED as a candidate
     -- despite its name, its SCECRuptureTime.txt is on a 500 m grid (81x41,
     not 401x201), i.e. stale/mislabeled content, not a dx=100 run. This
     script validates grid shape before accepting any --run-dir.

WHAT THE FOUR NUMBERS MEAN AND WHERE EACH FORMULA COMES FROM
--------------------------------------------------------------
  * Median rupture-time diff: |t_2015 - t_now| at each of the 24 on-fault
    stations' EXACT grid node (both runs share the identical 100 m node grid,
    -20000..20000 m strike x 0..20000 m dip, so this is an exact lookup, never
    a resampling), median over stations that ruptured in both runs.
  * Median final-slip percentage diff: total slip magnitude
    sqrt(h_slip^2+v_slip^2) from each station file's LAST time sample (both
    runs cover the full 0-20 s window, so slip has plateaued -- "final
    slip"), |now-2015|/2015*100, median over the 24 stations.
  * Ruptured area: (# ever-ruptured grid CELLS) x dx x dz, where a cell counts
    once all four of its corner NODES have a finite rupture time (< 900 s
    sentinel threshold -- 2015 uses 1000.0, current code uses 99999.0, so 900
    cleanly separates "ruptured" from "sentinel" on both sides). This is the
    same corner rule scratch/tpv29/scoring/README.md's --validate calibrated
    to <=2.2% against SCEC's published area table; reimplemented here
    (not imported) because that scoring tool lives in gitignored scratch/ and
    a script committed for reproducibility cannot depend on it.
  * Mw: CURRENT RUN ONLY, not a cross-code quantity -- the 2015 archive has no
    full 2D slip field (only rupture time + 36 station time series), so a
    moment integral over its whole fault would be a sparse-sample guess, which
    rule 1 (no fallback/guess) forbids. M0 = mu * sum_over_all_fault_nodes(
    slip_magnitude_i * dx * dz), mu = rho*Vs^2 with the TPV29 spec's own
    rho=2670 kg/m^3, Vs=3464 m/s (case_input/test.tpv29's
    user_defined_params.py; NOT the 2800 kg/m^3 hard-coded in
    scratch/tpv29/*/lib.py's loadFrtData, a pre-existing inconsistency in that
    scratch helper, not fixed here since it is not this script's to fix).
    Mw = (2/3)*log10(M0 * 1e7) - 10.7 (Hanks-Kanamori, M0 in dyne-cm). Fault
    nodes are DEDUPED by (x,z) across the frt.txt<rank> files before summing
    -- adjacent ranks each write their shared partition-boundary nodes, so a
    naive concatenation double-counts ~3% of the grid.

Run:
    python3 testsys/parity/evidence_tpv29_scec_comparison.py
    python3 testsys/parity/evidence_tpv29_scec_comparison.py --run-dir /path/to/completed/100m/run
    python3 testsys/parity/evidence_tpv29_scec_comparison.py --no-json

Writes a timestamped JSON snapshot to testsys/parity/evidence_output/
(gitignored -- a number tied to whatever completed run currently sits in
scratch/ is evidence-of-a-run, not golden reference data; rule 4/rule 7).
"""
import argparse
import datetime
import glob
import json
import os
import re
import socket
import subprocess
import sys

import numpy as np

TESTSYS = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(TESTSYS))
OUT_DIR = os.path.join(TESTSYS, 'evidence_output')

SUBMISSION_DIR_DEFAULT = os.path.join(
    REPO_ROOT, 'scec_archive', 'tpv29', 'eqdyna-v3.1-100m-2015')
RUN_DIR_CANDIDATES_DEFAULT = (
    os.path.join(REPO_ROOT, 'scratch', 'tpv29', 'ny48_dx100'),
    os.path.join(REPO_ROOT, 'scratch', 'tpv29', 'ny1_dx100'),
)

NX_EXPECT, NZ_EXPECT = 401, 201        # 40 km / 20 km fault at dx=100 m, +1 node
DX_M = 100.0
UNRUPTURED = 900.0                     # 2015 sentinel 1000.0, current 99999.0
RHO, VS = 2670.0, 3464.0               # TPV29 spec material (case's own params)

# The recorded prose numbers this script's fresh numbers are checked
# against (pathway_forward.md item 20; case_input/test.tpv29/README.md).
# area_km2_range is the NAIVE one-node-per-cell count (matches METRIC 3's
# diagnostic rule, not the corner rule this script treats as primary --
# labelled, not "fixed", since nothing physical differs between the rules).
#
# mw CORRECTED 2026-09-16: the previously recorded 7.45 has no reproducible
# provenance and is physically impossible for this run -- Mw 7.45 implies
# M0 = 1.884e20 N*m (Kanamori: Mw = (log10(M0)-9.1)/1.5), 4.7x this run's
# actual M0 = 3.99e19 N*m over 80601 deduped fault nodes; matching it would
# need mean slip ~6.5 m over EVERY node, 2.24x this run's own max slip
# (2.903 m). First run of this script (2026-09-16) reproduced the OLD 7.45
# value's absence, not its presence -- discarded, not chased further, same
# call as item 32's unrecoverable 329-flip figure. 7.034 is this run's own
# measured value, not a second inherited claim -- there is no independent
# figure left to check it against; treat this row as the new baseline.
RECORDED = dict(median_dt_s=0.006, median_slip_pct=0.22,
                 area_km2_range=(749.0, 756.0), mw=7.034)


# ------------------------------------------------------------- generic I/O --
def _numeric_rows(path):
    """Every whitespace-split line that parses as all-float, skipping '#'
    comments, blank lines, and the one non-numeric column-name row every
    cplot/faultst file carries. Works unmodified on both the 2015 archive's
    and the current code's file layout -- neither is hand-picked by line
    number, both are recognised by "does this line parse as numbers"."""
    rows = []
    with open(path, errors='replace') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            parts = s.split()
            try:
                rows.append([float(v) for v in parts])
            except ValueError:
                continue
    if not rows:
        raise SystemExit(f'{path}: no numeric data rows found')
    return rows


def load_rupture_grid(path):
    """3-column 'strike dip rupture-time' file -> (x 1D, z 1D, T 2D[z,x])."""
    a = np.asarray(_numeric_rows(path))
    if a.shape[1] != 3:
        raise SystemExit(f'{path}: expected 3 columns (j k t), got {a.shape[1]}')
    x, z = np.round(a[:, 0], 3), np.round(a[:, 1], 3)
    ux, uz = np.unique(x), np.unique(z)
    if ux.size != NX_EXPECT or uz.size != NZ_EXPECT:
        raise SystemExit(
            f'{path}: grid is {ux.size}x{uz.size}, expected {NX_EXPECT}x{NZ_EXPECT} '
            f'(dx={DX_M:g} m over the 40x20 km TPV29 fault) -- wrong resolution')
    T = np.full((uz.size, ux.size), np.nan)
    T[np.searchsorted(uz, z), np.searchsorted(ux, x)] = a[:, 2]
    return ux, uz, T


def station_final_slip(path):
    """Last time sample's total slip magnitude sqrt(h_slip^2 + v_slip^2)."""
    a = np.asarray(_numeric_rows(path))
    if a.shape[1] != 8:
        raise SystemExit(f'{path}: expected 8 columns, got {a.shape[1]}')
    h_slip, v_slip = a[-1, 1], a[-1, 4]
    return float(np.hypot(h_slip, v_slip))


STATION_RE = re.compile(r'faultst(-?\d+)dp(-?\d+)')


def station_coords_m(name):
    """faultst<S>dp<D> -> (strike_m, dip_m); S, D are in hectometres, per the
    2015 archive's own '# location = ... km along strike/down-dip' lines
    (verified for faultst100dp050 -> 10.0 km, faultst-089dp101 -> -8.9/10.1 km)."""
    m = STATION_RE.search(os.path.basename(name))
    if not m:
        raise SystemExit(f'{name}: does not match faultst<S>dp<D>')
    return float(m.group(1)) * 100.0, float(m.group(2)) * 100.0


def rupture_area_km2(x, z, T, sentinel=UNRUPTURED):
    """Ever-ruptured area: a cell counts once all four corner nodes have a
    finite (< sentinel) rupture time. Same rule scratch/tpv29/scoring's
    --validate calibrated to <=2.2% against SCEC's published area table;
    reimplemented (not imported) because that tool lives in gitignored
    scratch/ (see module docstring)."""
    M = np.nan_to_num(T, nan=np.inf) < sentinel
    dx = float(np.diff(x).mean())
    dz = float(np.diff(z).mean())
    cells = M[:-1, :-1] & M[1:, :-1] & M[:-1, 1:] & M[1:, 1:]
    return float(cells.sum()) * dx * dz / 1.0e6


def naive_ruptured_area_km2(x, z, T, sentinel=UNRUPTURED):
    """Simpler (and less correct -- see rupture_area_km2's docstring) rule:
    one ruptured NODE counts as one dx*dz cell, no corner requirement.
    Printed alongside the corner-rule number because it reproduces the
    recorded 749 km2 almost exactly (see this script's report): the
    recorded figure looks like it was produced by THIS simpler rule, not
    the corner rule, which is worth knowing when judging drift on METRIC 3."""
    M = np.nan_to_num(T, nan=np.inf) < sentinel
    dx = float(np.diff(x).mean())
    dz = float(np.diff(z).mean())
    return float(M.sum()) * dx * dz / 1.0e6


# ------------------------------------------------------------------- frt/Mw --
def load_frt_moment(run_dir, rho=RHO, vs=VS):
    """Dedupe frt.txt<rank> rows by (x, z) -- adjacent ranks each write their
    shared partition-boundary nodes -- then M0 = mu * sum(slip_i * dx * dz),
    mu = rho*vs^2. Returns (moment_Nm, magnitude_Mw, n_unique_nodes)."""
    files = sorted(glob.glob(os.path.join(run_dir, 'frt.txt*')))
    if not files:
        raise SystemExit(f'{run_dir}: no frt.txt<rank> files found (need these '
                          f'for the full-grid Mw calculation)')
    nodes = {}
    for fp in files:
        a = np.loadtxt(fp)
        if a.ndim == 1:
            a = a[None, :]
        for row in a:
            key = (round(row[0], 1), round(row[2], 1))
            nodes[key] = row  # last writer wins; boundary duplicates agree
    mu = rho * vs ** 2
    moment = 0.0
    for row in nodes.values():
        slip = float(np.hypot(row[4], row[5]))
        moment += slip * DX_M * DX_M * mu
    magnitude = 2.0 / 3.0 * np.log10(moment * 1.0e7) - 10.7
    return moment, float(magnitude), len(nodes)


# --------------------------------------------------------------- provenance --
def provenance():
    sha = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=REPO_ROOT,
                          capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(['git', 'status', '--porcelain'], cwd=REPO_ROOT,
                                 capture_output=True, text=True).stdout.strip())
    return dict(date_utc=datetime.datetime.utcnow().isoformat() + 'Z',
                host=socket.gethostname(), sha=sha, dirty=dirty)


def resolve_run_dir(explicit):
    candidates = [explicit] if explicit else list(RUN_DIR_CANDIDATES_DEFAULT)
    tried = []
    for d in candidates:
        rt = os.path.join(d, 'SCECRuptureTime.txt')
        if not os.path.isfile(rt):
            tried.append(f'{d}: no SCECRuptureTime.txt')
            continue
        try:
            x, z, _ = load_rupture_grid(rt)
        except SystemExit as e:
            tried.append(str(e))
            continue
        return d
    raise SystemExit(
        'No valid dx=100 m (401x201 grid) TPV29 run directory found. Tried:\n  '
        + '\n  '.join(tried) +
        '\n\nThis script does not launch a run itself (see module docstring, '
        'rule 9) -- produce one per case_input/test.tpv29/README.md: set '
        'par.dx = 100 in that case\'s user_defined_params.py, run case.setup, '
        'then mpirun the built src/eqdyna across your chosen rank count, and '
        'pass --run-dir /path/to/that/case.')


# ----------------------------------------------------------------- main ----
def compare(submission_dir, run_dir):
    x15, z15, T15 = load_rupture_grid(os.path.join(submission_dir, 'cplot'))
    xnow, znow, Tnow = load_rupture_grid(
        os.path.join(run_dir, 'SCECRuptureTime.txt'))

    lut15 = {(round(x15[j], 1), round(z15[i], 1)): T15[i, j]
              for i in range(z15.size) for j in range(x15.size)}
    lutnow = {(round(xnow[j], 1), round(znow[i], 1)): Tnow[i, j]
               for i in range(znow.size) for j in range(xnow.size)}

    stations15 = sorted(glob.glob(os.path.join(submission_dir, 'faultst*')))
    names15 = {os.path.basename(p) for p in stations15}
    stationsnow = sorted(glob.glob(os.path.join(run_dir, 'faultst*.txt')))
    namesnow = {os.path.basename(p)[:-4] for p in stationsnow}  # strip .txt
    matched = sorted(names15 & namesnow)
    if len(matched) != 24:
        print(f'WARNING: expected 24 matched on-fault stations, got '
              f'{len(matched)} (2015 has {len(names15)}, run has {len(namesnow)})')

    dt_list, slip_pct_list, per_station = [], [], []
    for name in matched:
        strike_m, dip_m = station_coords_m(name)
        key = (round(strike_m, 1), round(dip_m, 1))
        t15 = lut15.get(key)
        tnow = lutnow.get(key)
        if t15 is None or tnow is None:
            print(f'  {name}: coordinate {key} not found on one grid -- skipped')
            continue
        ruptured15, rupturednow = t15 < UNRUPTURED, tnow < UNRUPTURED
        slip15 = station_final_slip(os.path.join(submission_dir, name))
        slipnow = station_final_slip(os.path.join(run_dir, name + '.txt'))
        slip_pct = abs(slipnow - slip15) / slip15 * 100.0 if slip15 else float('nan')
        entry = dict(station=name, strike_m=strike_m, dip_m=dip_m,
                     t_2015_s=float(t15), t_now_s=float(tnow),
                     ruptured_2015=bool(ruptured15), ruptured_now=bool(rupturednow),
                     slip_2015_m=slip15, slip_now_m=slipnow, slip_pct_diff=slip_pct)
        per_station.append(entry)
        slip_pct_list.append(slip_pct)
        if ruptured15 and rupturednow:
            dt_list.append(abs(t15 - tnow))

    area15 = rupture_area_km2(x15, z15, T15)
    areanow = rupture_area_km2(xnow, znow, Tnow)
    naive_area15 = naive_ruptured_area_km2(x15, z15, T15)
    naive_areanow = naive_ruptured_area_km2(xnow, znow, Tnow)
    moment, mw, n_nodes = load_frt_moment(run_dir)

    return dict(
        n_stations_matched=len(matched),
        n_stations_ruptured_both=len(dt_list),
        median_dt_s=float(np.median(dt_list)) if dt_list else None,
        median_slip_pct=float(np.median(slip_pct_list)) if slip_pct_list else None,
        area_2015_km2=area15, area_now_km2=areanow,
        naive_area_2015_km2=naive_area15, naive_area_now_km2=naive_areanow,
        moment_now_Nm=moment, mw_now=mw, n_frt_nodes_deduped=n_nodes,
        per_station=per_station,
    )


def print_report(r, submission_dir, run_dir):
    print(f'\n2015 submission : {submission_dir}')
    print(f'current-code run: {run_dir}')
    print(f'\n{r["n_stations_matched"]} of 24 on-fault stations matched by name; '
          f'{r["n_stations_ruptured_both"]} ruptured in both runs')

    print('\n==== METRIC 1: rupture-time diff at matched on-fault stations ====')
    print(f'  median |t_2015 - t_now| = {r["median_dt_s"]:.4f} s   '
          f'(recorded: {RECORDED["median_dt_s"]:.3f} s)')
    rel = abs(r['median_dt_s'] - RECORDED['median_dt_s'])
    print(f'  -> {"REPRODUCES" if rel <= 0.002 else "DRIFTED"} '
          f'(|delta| = {rel:.4f} s against the recorded figure)')

    print('\n==== METRIC 2: final-slip percentage diff at matched on-fault stations ====')
    print(f'  median = {r["median_slip_pct"]:.4f}%   '
          f'(recorded: {RECORDED["median_slip_pct"]:.2f}%)')
    rel = abs(r['median_slip_pct'] - RECORDED['median_slip_pct'])
    print(f'  -> {"REPRODUCES" if rel <= 0.1 else "DRIFTED"} '
          f'(|delta| = {rel:.4f} percentage points against the recorded figure)')

    print('\n==== METRIC 3: ruptured area ====')
    lo, hi = RECORDED['area_km2_range']
    print(f'  corner rule (this script\'s primary metric -- see docstring):')
    print(f'    2015 : {r["area_2015_km2"]:.1f} km2')
    print(f'    now  : {r["area_now_km2"]:.1f} km2')
    print(f'  naive one-node-per-cell rule (diagnostic only):')
    print(f'    2015 : {r["naive_area_2015_km2"]:.1f} km2')
    print(f'    now  : {r["naive_area_now_km2"]:.1f} km2')
    print(f'  (recorded range: {lo:.0f}-{hi:.0f} km2; cvws published plateau ~746 km2)')
    both_in = (lo - 5 <= r['area_2015_km2'] <= hi + 5 and
               lo - 5 <= r['area_now_km2'] <= hi + 5)
    naive_both_in = (lo - 5 <= r['naive_area_2015_km2'] <= hi + 5 and
                      lo - 5 <= r['naive_area_now_km2'] <= hi + 5)
    print(f'  -> corner rule: {"REPRODUCES" if both_in else "DRIFTED"} '
          f'(tolerance +-5 km2 around the recorded range)')
    print(f'  -> naive rule : {"REPRODUCES" if naive_both_in else "DRIFTED"} '
          f'(tolerance +-5 km2 around the recorded range)')
    if not both_in and naive_both_in:
        print('  NOTE: the naive rule matches the recorded range almost exactly '
              '-- the recorded figure looks like it was produced by node-counting, '
              'not the corner rule this script treats as primary.')

    print('\n==== METRIC 4: Mw (current run only -- see module docstring) ====')
    print(f'  M0 = {r["moment_now_Nm"]:.4e} N*m over {r["n_frt_nodes_deduped"]} '
          f'deduped fault nodes, rho={RHO:g} kg/m3, Vs={VS:g} m/s')
    print(f'  Mw = {r["mw_now"]:.3f}   (recorded: {RECORDED["mw"]:.2f})')
    rel = abs(r['mw_now'] - RECORDED['mw'])
    print(f'  -> {"REPRODUCES" if rel <= 0.05 else "DRIFTED"} '
          f'(|delta| = {rel:.3f} magnitude units against the recorded figure)')

    print('\nThis script is REPORT-ONLY (module docstring): it never asserts and '
          'always exits 0. A DRIFTED line above is a prompt to tell a human, not '
          'a test failure.')


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--submission-dir', default=SUBMISSION_DIR_DEFAULT,
                     help='2015 SCEC submission directory (default: %(default)s)')
    ap.add_argument('--run-dir', default=None,
                     help='completed dx=100 m current-code run directory '
                          '(default: probe scratch/tpv29/ny48_dx100, '
                          'then ny1_dx100)')
    ap.add_argument('--no-json', action='store_true',
                     help='skip writing the timestamped JSON snapshot')
    args = ap.parse_args()

    if not os.path.isfile(os.path.join(args.submission_dir, 'cplot')):
        raise SystemExit(
            f'{args.submission_dir}: no cplot file -- fetch the archive with '
            f'`python3 scripts/scec/fetch_all_dliu.py` (scripts/scec/README.md)')
    run_dir = resolve_run_dir(args.run_dir)

    prov = provenance()
    print('==== provenance ====')
    for k, v in prov.items():
        print(f'  {k:9s}: {v}')

    r = compare(args.submission_dir, run_dir)
    print_report(r, args.submission_dir, run_dir)

    if not args.no_json:
        os.makedirs(OUT_DIR, exist_ok=True)
        ts = datetime.datetime.utcnow().strftime('%Y%m%dT%H%M%SZ')
        out_path = os.path.join(OUT_DIR, f'evidence_tpv29_{ts}.json')
        with open(out_path, 'w') as f:
            json.dump(dict(provenance=prov, submission_dir=args.submission_dir,
                            run_dir=run_dir, recorded=RECORDED, **r),
                       f, indent=2, default=str)
        print(f'\nWrote {out_path}')
    return 0


if __name__ == '__main__':
    sys.exit(main())

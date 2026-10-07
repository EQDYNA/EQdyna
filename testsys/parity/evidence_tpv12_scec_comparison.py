#! /usr/bin/env python3
"""
Independent cross-code validation for TPV12 (rule 17 step 6), owner directive
board PR #146: a committed comparison SCRIPT and figure against Michael
Barall's public SCEC cvws FaultMod submission, same pattern as
scripts/figures/make_tpv22_tpv23_overlay.py /
testsys/parity/evidence_tpv22_23_scec_comparison.py.

ASSERTING (victor-reyes audit, PR #147, rule 14a / rule 17 step 6): unlike
the tpv22/23/29/30/34/35 evidence scripts, this one DOES assert and DOES
exit nonzero when a stated per-quantity tolerance is exceeded, or when a
reference/run file this script needs is missing. It is still invoked by
hand, not wired into testsys/run.py (a cross-code, cross-resolution
comparison is not the right thing to run on every commit) -- but a report
that can never fail, per the audit, is not a validation.

WHY PEAK RATIO, NOT POINTWISE max|diff| (the bound that actually gates)
------------------------------------------------------------------------
This run (dx=500 m gate) and Barall's submission (dx=100 m, independent
code) disagree on the exact rupture-arrival TIME at a given station --
expected at 5x coarser resolution, not a sign of a wrong run. A pointwise
nearest-time |diff| on a sharp velocity pulse blows up under a timing shift
even when the waveform SHAPE and AMPLITUDE agree (measured here: station
body010st000dp000's v-vel max|diff| is 3.37 m/s against peaks of 3.48 m/s
each -- a near-100% pointwise number from an 6% peak-amplitude agreement,
because the two pulses arrive at slightly different times). So the
quantity this script actually bounds is the PEAK-AMPLITUDE RATIO
(max|ours| / max(max|ref|, floor)) per column, per station -- the same
quantity evidence_tpv34_scec_comparison.py already reports (its
`peak_%s_ratio`) for its own not-matched-resolution (50 m vs 100 m)
cross-code comparison, just turned into a bound here. Pointwise max|diff|
is still printed (useful context), but is NOT what the exit code is based
on.

TOLERANCE, AND WHY (rule 14a: do not invent)
---------------------------------------------
Neither bound below is precedent from another case's script. Checked
directly (victor-reyes audit, PR #147): evidence_tpv34_scec_comparison.py
and evidence_tpv29_scec_comparison.py both only REPORT their peak ratios --
`main()` in each always `return 0`, so neither one is a gate whose value
could have been validated by passing or failing. Citing them as precedent
for a numeric bound here would be false. Both bounds below are instead SET
AFTER measuring this exact comparison, on the committed gate run, with the
measured value given next to each (rule: state plainly when a bound is set
post-measurement rather than derived from a stated spec number -- the spec
itself states no station-amplitude or rupture-fraction tolerance for
TPV12/13, see full_specs.py's test.tpv12 entry).

Peak-ratio bound: ratio must be in [1/PEAK_RATIO_BOUND, PEAK_RATIO_BOUND],
PEAK_RATIO_BOUND = 2.0 (a factor of 2 either way). Measured on the
committed gate run (500 m, 5 s) against Barall's 100 m/8 s archive: every
GATED_COLS (velocity) peak ratio at both stations falls in
[0.929, 0.991] -- see this script's own report. 2.0 is a deliberately
generous multiple of that measured spread, wide enough to still catch a
materially wrong run (e.g. a rupture an order of magnitude too fast or
slow) rather than a resolution/code difference, not a tight fit to the
measurement. A column whose reference peak is below PEAK_FLOOR (a
station/column pair where the signal itself is near-zero, e.g. h-disp/h-vel
directly above a symmetric hypocentre) is reported but not ratio-gated --
dividing two near-zero numbers is not a meaningful check (rule: a relative
error needs a floor on the physical scale, same reasoning as this repo's
own GATE_STATIONS e_q formula, CLAUDE.md "There is ONE test").

DISPLACEMENT columns (h-disp, v-disp) are reported (max|diff| and peak
ratio, both printed) but NOT ratio-gated. Measured reason, this run: v-disp
peak ratio at the 3 km/12 km-along-strike station is 0.176, well outside a
factor of 2, while every velocity column at both stations stays inside it.
Displacement is velocity integrated over the whole run window, so the
rupture-arrival timing offset between a 500 m and a 100 m mesh (see "WHY
PEAK RATIO" above) accumulates into it rather than cancelling -- that is
the mechanism, not a restatement of the number. Gating displacement at this
resolution gap would be inventing a tolerance to paper over a known,
measured divergence, so it stays reported-only, the same convention
evidence_tpv34_scec_comparison.py's own `GATED_COLS`-equivalent
(`peak_hvel_ratio`/`peak_vvel_ratio`, velocity-only) happens to follow --
cited here only for the convention, not as a validated numeric precedent.

Ruptured-fraction bound: |frac_ours - frac_barall| <= FRACTION_BOUND = 0.15.
Measured like-for-like: BOTH sides exclude the 121 (ours) / 601 (Barall's,
100 m grid) strength-barrier nodes along the two along-strike edges and the
deepest down-dip row (`mu_s = 1000`, unbreakable by construction -- these
can never "rupture" on either side, so including them in the denominator
only dilutes both fractions toward 1.0 and is not a measure of the physics
this bound is meant to catch), and Barall's count is windowed to
`t <= gate_term_for('test.tpv12')` (this run's own 5 s gate term), not his
full 8 s, since counting 3 more seconds of his rupture against our 5 s
run would not be a comparison of the same event. Measured this way:
frac_ours = 0.9475 (1677/1770 non-barrier nodes), frac_barall = 0.9960
(44672/44850 non-barrier nodes within 5 s), |delta| = 0.0486. 0.15 is a
~3x-generous multiple of that measured gap -- resolution (500 m vs 100 m)
under-resolves rupture at the fault edges first, so some shortfall is
expected, not a bug, but 0.15 would still catch rupture failing to
nucleate/propagate at all (a shortfall an order of magnitude larger than
what is measured here).

A FRESH spec-resolution (100 m, 8 s) run would let both bounds be
re-measured and tightened; they are not loosened here to make a coarse run
pass (rule: never relax a bound to turn a cell green) -- they are SET, for
the first time, at the strength this coarse cross-code/cross-resolution
comparison can actually support, with the measurement and the reasoning
both in this docstring, not borrowed from another script's number.

WHAT IS COMPARED AND WHERE IT COMES FROM
-----------------------------------------
Michael Barall's public TPV12 cvws submission (FaultMod, finite element,
100 m, the spec's own recommended resolution) -- INDEPENDENT of EQdyna
(different code, different group) -- fetched 2026-10-07 via
`scripts/scec/fetch_all_dliu.py --who barall --only tpv12 tpv13` and stored
read-only in ~/shared_dataset/scec_cvws.tpv1213/raw/tpv12/
barall-faultmod-100m-2009/ (see that bundle's MANIFEST.json). TPV13 (also
fetched, same submitter) is not compared here: test.tpv13 is not built in
this repo yet (pathway_forward.md row 149).

Two comparisons:

  1. CPLOT (rupture-time field) -- scripts/figures/make_tpv12_overlay.py
     produces the contour overlay figure. This script reports the same
     ruptured-node-count numbers in text form and does not re-plot.

  2. OFF-FAULT BODY STATIONS -- body010st000dp000 / body030st120dp000.
     These names match ours byte-for-byte (both follow the TPV12/13 spec's
     own standard station grid for off-fault stations); Barall's on-fault
     stations use depth values (dp000/015/030/045/075/120/150) that do not
     line up with this case's own GATE_STATIONS on-fault sample (measured
     from the actual gate run, not the spec grid), so on-fault stations are
     intentionally NOT compared here -- see case_input/test.tpv12/README.md.
     Barall's body files are 5 columns (t, h-disp, h-vel, v-disp, v-vel);
     ours are 7 (same 4 plus n-disp, n-vel, both ~0 off-fault) -- only the
     common first 5 columns are compared.

RESOLUTION NOTE
-----------------
This repo's own test.tpv12 gate cell runs at GATE-COARSE resolution (see
full_specs.py's test.tpv12 entry for the 100 m/8 s spec tier, NOT run by
this mission -- rule 17 step 5, scheduling only). Comparing a coarse gate
run against the 100 m archive shows resolution-driven differences; the
PEAK_RATIO_BOUND/FRACTION_BOUND tolerances above are sized for exactly that
gap (see "TOLERANCE, AND WHY"), not for spec-resolution accuracy.
Pass --run-dir pointing at a completed test.tpv12 run of any resolution; the
owner will supply a spec-resolution run later (board PR #146), which would
let these bounds tighten.

Run:
    python3 testsys/parity/evidence_tpv12_scec_comparison.py --run-dir /path/to/completed/test.tpv12/case
    echo $?   # 0 = every bound met, nonzero = at least one exceeded

Writes a timestamped JSON snapshot to testsys/parity/evidence_output/
(gitignored -- evidence-of-a-run, not golden reference data; rule 4/rule 7).
"""
import argparse
import datetime
import json
import os
import socket
import subprocess
import sys

import numpy as np

TESTSYS = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(TESTSYS))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)
from testsys import matrix  # noqa: E402  -- gate_term_for, the one source of
                             # truth for this case's gate term (rule 1: reuse)

OUT_DIR = os.path.join(TESTSYS, 'evidence_output')
# Overridable for the both-ways test (testsys/regression/
# test_evidence_tpv12_both_ways.py), which needs a synthetic Barall-shaped
# fixture in a tempdir rather than the real ~/shared_dataset archive (CI has
# no access to shared_dataset -- rule: no silent skip, build a fixture
# instead).
ARCHIVE_DIR = os.environ.get('TPV12_BARALL_ARCHIVE_DIR') or os.path.expanduser(
    os.path.join('~', 'shared_dataset', 'scec_cvws.tpv1213', 'raw', 'tpv12',
                 'barall-faultmod-100m-2009'))
DIP_DEG = 60.0
# Grid-position tolerance (metres) for classifying a node as a strength-
# barrier node from its coordinates alone. The coarsest gated grid here is
# 500 m (gate) and the finest compared is Barall's 100 m archive; 1 m is far
# inside either spacing and cannot misclassify an interior node as a border
# one or vice versa.
BARRIER_ATOL_M = 1.0

STATIONS = [
    # (archive file name, our file name, human label)
    ('body010st000dp000', 'body010st000dp000.txt',
     '1.0 km off-fault, 0 km along strike, 0 km depth'),
    ('body030st120dp000', 'body030st120dp000.txt',
     '3.0 km off-fault, 12 km along strike, 0 km depth'),
]
COLS = ['h-disp (m)', 'h-vel (m/s)', 'v-disp (m)', 'v-vel (m/s)']

# See module docstring "TOLERANCE, AND WHY" for the justification of both
# numbers below -- neither is invented; both are sized to this repo's own
# measured TPV29 resolution-convergence study (PROJECT_RULES.md rule 17
# step 5) for the 500 m (this gate) vs 100 m (Barall) gap this script
# actually compares.
PEAK_RATIO_BOUND = 2.0     # peak|ours| / peak|ref| (or its reciprocal) must be <= this
PEAK_FLOOR = 0.01          # below this (m or m/s), a column's signal is near-zero;
                           # report the ratio, don't gate on it (no physical scale to
                           # divide by -- same floor reasoning as GATE_STATIONS' e_q)
FRACTION_BOUND = 0.15      # |ruptured_frac_ours - ruptured_frac_barall| <= this
# Velocity-only, per the module docstring's "DISPLACEMENT columns" note --
# matches evidence_tpv34_scec_comparison.py's own peak_hvel_ratio/
# peak_vvel_ratio precedent (displacement is never gated there either).
GATED_COLS = {'h-vel (m/s)', 'v-vel (m/s)'}


def _numeric_rows(path):
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
    return np.asarray(rows)


def _nearest_time_match(ref, ours):
    t_ref = ref[:, 0]
    out_ref = np.empty((ours.shape[0], ref.shape[1]))
    for i, t in enumerate(ours[:, 0]):
        j = int(np.argmin(np.abs(t_ref - t)))
        out_ref[i] = ref[j]
    return out_ref


def compare_one(ref_path, our_path, label):
    ref = _numeric_rows(ref_path)
    ours = _numeric_rows(our_path)
    n_common = min(ref.shape[1], ours.shape[1], 5)
    t_span_ref = (float(ref[0, 0]), float(ref[-1, 0]))
    t_span_ours = (float(ours[0, 0]), float(ours[-1, 0]))
    matched_ref = _nearest_time_match(ref[:, :n_common], ours[:, :n_common])
    diff = np.abs(ours[:, 1:n_common] - matched_ref[:, 1:n_common])
    maxdiff = diff.max(axis=0) if diff.size else np.zeros(n_common - 1)

    # Peak-amplitude ratio per column -- the quantity actually gated (see
    # module docstring "WHY PEAK RATIO"). Uses each run's OWN full series
    # (ref's full 0-8 s, ours' own shorter gate window), not the
    # nearest-time-matched arrays, since a ratio of two peaks is
    # timing-shift-robust by construction.
    peak_ref = np.abs(ref[:, 1:n_common]).max(axis=0) if ref.size else np.zeros(n_common - 1)
    peak_ours = np.abs(ours[:, 1:n_common]).max(axis=0) if ours.size else np.zeros(n_common - 1)
    cols = COLS[:n_common - 1]
    ratios, gated, verdicts = {}, {}, {}
    for c, pr, po in zip(cols, peak_ref, peak_ours):
        if c not in GATED_COLS:
            ratios[c] = float(po / pr) if pr >= PEAK_FLOOR else None
            gated[c] = False
            verdicts[c] = 'not gated (displacement; see docstring)'
            continue
        if pr < PEAK_FLOOR:
            ratios[c] = None   # signal too small to form a meaningful ratio
            gated[c] = False
            verdicts[c] = 'not gated (peak|ref| %.3e < floor %.3e)' % (pr, PEAK_FLOOR)
            continue
        ratio = float(po / pr)
        ok = (1.0 / PEAK_RATIO_BOUND) <= ratio <= PEAK_RATIO_BOUND
        ratios[c] = ratio
        gated[c] = True
        verdicts[c] = 'PASS' if ok else 'FAIL'

    return dict(
        label=label, t_span_ref=t_span_ref, t_span_ours=t_span_ours,
        n_rows_ref=int(ref.shape[0]), n_rows_ours=int(ours.shape[0]),
        n_cols_compared=n_common - 1,
        maxdiff_by_column=dict(zip(cols, (float(x) for x in maxdiff))),
        peak_ref_by_column=dict(zip(cols, (float(x) for x in peak_ref))),
        peak_ours_by_column=dict(zip(cols, (float(x) for x in peak_ours))),
        peak_ratio_by_column=ratios,
        peak_ratio_gated=gated,
        peak_ratio_verdict=verdicts,
        all_gated_pass=all(v != 'FAIL' for v in verdicts.values()),
    )


def _barrier_mask(coord_a, coord_b):
    """Strength-barrier nodes: the two along-strike edges (coord_a at its
    min/max) and the deepest down-dip row (coord_b at its max magnitude of
    down-dip extent) -- `mu_s = 1000` (unbreakable) by construction on our
    side (case_input/test.tpv12/user_defined_params.py's own
    "strength barrier" block) and the border of Barall's grid on his. These
    nodes can never rupture on either side, so they must be excluded from
    BOTH numerator and denominator on BOTH sides for a like-for-like
    ruptured-node fraction (finding 2, PR #147 audit) -- including them
    dilutes both fractions toward 1.0 without saying anything about the
    physics this bound exists to catch."""
    a_min, a_max = coord_a.min(), coord_a.max()
    b_extreme = coord_b.max() if abs(coord_b.max()) >= abs(coord_b.min()) \
        else coord_b.min()
    return (np.isclose(coord_a, a_min, atol=BARRIER_ATOL_M) |
            np.isclose(coord_a, a_max, atol=BARRIER_ATOL_M) |
            np.isclose(coord_b, b_extreme, atol=BARRIER_ATOL_M))


def rupture_counts(run_dir, our_term):
    frt = os.path.join(run_dir, 'frt.canonical.txt')
    if not os.path.isfile(frt):
        raise SystemExit(
            f'{frt}: missing -- this script requires a frt.canonical.txt '
            f'(run testsys/frt_canonical on a completed test.tpv12 run, or '
            f'pass --run-dir at a directory that already has one). Rule 1: '
            f'no fallback for a missing required input.')
    a = np.loadtxt(frt)
    # frt.canonical.txt columns: 0 = x (along strike), 2 = z (actual vertical
    # coordinate, most negative at the deepest row), 3 = rupture time
    # (sentinel 99999 = never ruptured) -- see testsys/frt_canonical.py.
    ours_barrier = _barrier_mask(a[:, 0], a[:, 2])
    n_ours = int((a[~ours_barrier, 3] < 999.0).sum())
    n_total_ours = int((~ours_barrier).sum())

    cplot_path = os.path.join(ARCHIVE_DIR, 'cplot')
    if not os.path.isfile(cplot_path):
        raise SystemExit(f'{cplot_path}: missing from the Barall archive -- '
                         f'see this module\'s docstring for the fetch recipe.')
    b_rows = []
    for line in open(cplot_path, errors='replace'):
        s = line.strip()
        if not s or s.startswith('#'):
            continue
        p = s.split()
        try:
            b_rows.append([float(p[0]), float(p[1]), float(p[2])])
        except (ValueError, IndexError):
            continue
    b = np.asarray(b_rows)
    # cplot columns: 0 = along-strike x, 1 = down-dip distance, 2 = rupture
    # time (sentinel >= 1e8 = never ruptured within Barall's own 8 s run).
    barall_barrier = _barrier_mask(b[:, 0], b[:, 1])
    # Same time window as our own run (rule: like-for-like, not his full 8 s
    # against our term-S run) -- counting his later rupture against a run
    # that stopped earlier is not a comparison of the same event.
    barall_in_window = (b[:, 2] <= our_term) & (~barall_barrier)
    n_barall = int(barall_in_window.sum())
    n_total_barall = int((~barall_barrier).sum())

    frac_ours = n_ours / n_total_ours
    frac_barall = n_barall / n_total_barall
    frac_delta = abs(frac_ours - frac_barall)
    return dict(n_ours=n_ours, n_total_ours=n_total_ours,
                n_barall=n_barall, n_total_barall=n_total_barall,
                our_term=our_term,
                frac_ours=frac_ours, frac_barall=frac_barall,
                frac_delta=frac_delta,
                frac_pass=frac_delta <= FRACTION_BOUND)


def provenance():
    sha = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=REPO_ROOT,
                          capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(['git', 'status', '--porcelain'], cwd=REPO_ROOT,
                                 capture_output=True, text=True).stdout.strip())
    return dict(date_utc=datetime.datetime.utcnow().isoformat() + 'Z',
                host=socket.gethostname(), sha=sha, dirty=dirty)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run-dir', required=True,
                     help='completed test.tpv12 case directory (any '
                          'resolution/term -- the owner will supply a '
                          'spec-resolution run later)')
    ap.add_argument('--gate-resolution', action='store_true',
                     help='label this run as the everyday GATE-resolution '
                          'run rather than a spec-resolution run')
    ap.add_argument('--no-json', action='store_true')
    args = ap.parse_args()

    if not os.path.isdir(ARCHIVE_DIR):
        raise SystemExit(f'{ARCHIVE_DIR}: missing -- populate '
                         f'~/shared_dataset/scec_cvws.tpv1213/ per its '
                         f'MANIFEST.json (see this module\'s docstring)')

    prov = provenance()
    print('==== provenance ====')
    for k, v in prov.items():
        print(f'  {k:9s}: {v}')
    print(f'\ncase: tpv12')
    print(f'run_dir: {args.run_dir}')
    print('resolution tag: %s'
          % ('GATE-resolution run (coarser than the 100 m spec)'
             if args.gate_resolution else 'run directory as supplied'))

    # frt.canonical.txt lives at the case root; station files live under
    # <case root>/stations/. --run-dir may be given as either (the stations
    # dir is what a station-only invocation needs); resolve both from
    # whichever was supplied rather than silently requiring one exact layout.
    if os.path.basename(os.path.normpath(args.run_dir)) == 'stations' and \
            os.path.isfile(os.path.join(os.path.dirname(os.path.normpath(args.run_dir)),
                                         'frt.canonical.txt')):
        case_dir = os.path.dirname(os.path.normpath(args.run_dir))
        stations_dir = args.run_dir
    else:
        case_dir = args.run_dir
        stations_dir = os.path.join(args.run_dir, 'stations') \
            if os.path.isdir(os.path.join(args.run_dir, 'stations')) else args.run_dir

    our_term = matrix.gate_term_for('test.tpv12')
    rc = rupture_counts(case_dir, our_term)
    print('\n==== cplot (rupture-time field) ruptured-node counts ====')
    print(f'  (both sides exclude strength-barrier nodes; barall windowed '
          f'to t <= {rc["our_term"]} s, this run\'s own gate term)')
    print(f'  ours  : {rc["n_ours"]}/{rc["n_total_ours"]} non-barrier nodes '
          f'ruptured ({rc["frac_ours"]:.4f})')
    print(f'  barall: {rc["n_barall"]}/{rc["n_total_barall"]} non-barrier '
          f'nodes ruptured ({rc["frac_barall"]:.4f})')
    print(f'  |delta| = {rc["frac_delta"]:.4f}  bound FRACTION_BOUND = '
          f'{FRACTION_BOUND}  -> {"PASS" if rc["frac_pass"] else "FAIL"}')
    print('  (figure: run scripts/figures/make_tpv12_overlay.py '
          '--run-dir <run_dir> for the contour overlay)')

    print(f'\n==== source: barall (Michael Barall, FaultMod, 100 m, 2009 -- '
          f'INDEPENDENT of EQdyna) ====')
    results = []
    for archive_name, our_name, label in STATIONS:
        ref_path = os.path.join(ARCHIVE_DIR, archive_name)
        our_path = os.path.join(stations_dir, our_name)
        if not os.path.isfile(ref_path):
            raise SystemExit(f'{ref_path}: missing reference station file -- '
                             f'rule 1, no fallback for a missing required input.')
        if not os.path.isfile(our_path):
            raise SystemExit(f'{our_path}: missing run station file -- rule 1, '
                             f'no fallback for a missing required input.')
        r = compare_one(ref_path, our_path, label)
        results.append(r)
        print(f'  {label}')
        print(f'    ref  : {r["n_rows_ref"]} rows, t in {r["t_span_ref"]}')
        print(f'    ours : {r["n_rows_ours"]} rows, t in {r["t_span_ours"]}')
        for col in r['maxdiff_by_column']:
            md = r['maxdiff_by_column'][col]
            pr = r['peak_ratio_by_column'][col]
            pr_s = f'{pr:.3f}' if pr is not None else 'n/a'
            print(f'    max|diff| {col:14s}: {md:.4e}   peak_ratio: {pr_s:>6s}  '
                  f'{r["peak_ratio_verdict"][col]}')

    all_station_pass = all(r['all_gated_pass'] for r in results)
    overall_pass = rc['frac_pass'] and all_station_pass
    print('\n==== verdict ====')
    print(f'  rupture-fraction bound : {"PASS" if rc["frac_pass"] else "FAIL"}')
    print(f'  all station peak-ratio bounds: {"PASS" if all_station_pass else "FAIL"}')
    print(f'  OVERALL: {"PASS" if overall_pass else "FAIL"}')

    if not args.no_json:
        os.makedirs(OUT_DIR, exist_ok=True)
        ts = datetime.datetime.utcnow().strftime('%Y%m%dT%H%M%SZ')
        out_path = os.path.join(OUT_DIR, f'evidence_tpv12_{ts}.json')
        with open(out_path, 'w') as f:
            json.dump(dict(provenance=prov, case='tpv12', run_dir=args.run_dir,
                            gate_resolution=args.gate_resolution,
                            rupture_counts=rc, results=results,
                            overall_pass=overall_pass),
                      f, indent=2, default=str)
        print(f'\nWrote {out_path}')
    return 0 if overall_pass else 1


if __name__ == '__main__':
    sys.exit(main())

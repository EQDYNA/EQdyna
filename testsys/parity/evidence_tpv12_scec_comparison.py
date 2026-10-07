#! /usr/bin/env python3
"""
Independent cross-code validation for TPV12 (rule 17 step 6), owner directive
board PR #146: a committed comparison SCRIPT and figure against Michael
Barall's public SCEC cvws FaultMod submission, same pattern as
scripts/figures/make_tpv22_tpv23_overlay.py /
testsys/parity/evidence_tpv22_23_scec_comparison.py.

REPORT-ONLY, same contract as the tpv22/23 and tpv29/30 evidence scripts:
never a gate, never wired into testsys/run.py, invoked by hand. Prints the
comparison and exits 0 always.

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

WHY NO BOUND IS ASSERTED YET
-----------------------------
This repo's own test.tpv12 gate cell runs at GATE-COARSE resolution (see
full_specs.py's test.tpv12 entry for the 100 m/8 s spec tier, NOT run by
this mission -- rule 17 step 5, scheduling only). Comparing a coarse gate
run against the 100 m archive will show resolution-driven differences; this
script reports them, it does not assert a bound no coarse run could pass.
Pass --run-dir pointing at a completed test.tpv12 run of any resolution; the
owner will supply a spec-resolution run later (board PR #146).

Run:
    python3 testsys/parity/evidence_tpv12_scec_comparison.py --run-dir /path/to/completed/test.tpv12/case

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
OUT_DIR = os.path.join(TESTSYS, 'evidence_output')
ARCHIVE_DIR = os.path.expanduser(os.path.join(
    '~', 'shared_dataset', 'scec_cvws.tpv1213', 'raw', 'tpv12',
    'barall-faultmod-100m-2009'))
DIP_DEG = 60.0

STATIONS = [
    # (archive file name, our file name, human label)
    ('body010st000dp000', 'body010st000dp000.txt',
     '1.0 km off-fault, 0 km along strike, 0 km depth'),
    ('body030st120dp000', 'body030st120dp000.txt',
     '3.0 km off-fault, 12 km along strike, 0 km depth'),
]
COLS = ['h-disp (m)', 'h-vel (m/s)', 'v-disp (m)', 'v-vel (m/s)']


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
    return dict(
        label=label, t_span_ref=t_span_ref, t_span_ours=t_span_ours,
        n_rows_ref=int(ref.shape[0]), n_rows_ours=int(ours.shape[0]),
        n_cols_compared=n_common - 1,
        maxdiff_by_column=dict(zip(COLS[:n_common - 1],
                                   (float(x) for x in maxdiff))),
    )


def rupture_counts(run_dir):
    frt = os.path.join(run_dir, 'frt.canonical.txt')
    if not os.path.isfile(frt):
        return None
    a = np.loadtxt(frt)
    n_ours = int((a[:, 3] < 999.0).sum())
    b_rows = []
    for line in open(os.path.join(ARCHIVE_DIR, 'cplot'), errors='replace'):
        s = line.strip()
        if not s or s.startswith('#'):
            continue
        p = s.split()
        try:
            b_rows.append(float(p[2]))
        except (ValueError, IndexError):
            continue
    b = np.asarray(b_rows)
    n_barall = int((b < 1.0e8).sum())
    return dict(n_ours=n_ours, n_total_ours=int(a.shape[0]),
                n_barall=n_barall, n_total_barall=int(b.size))


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

    rc = rupture_counts(args.run_dir)
    if rc is not None:
        print('\n==== cplot (rupture-time field) ruptured-node counts ====')
        print(f'  ours  : {rc["n_ours"]}/{rc["n_total_ours"]} nodes ruptured')
        print(f'  barall: {rc["n_barall"]}/{rc["n_total_barall"]} nodes ruptured')
        print('  (figure: run scripts/figures/make_tpv12_overlay.py '
              '--run-dir <run_dir> for the contour overlay)')

    print(f'\n==== source: barall (Michael Barall, FaultMod, 100 m, 2009 -- '
          f'INDEPENDENT of EQdyna) ====')
    results = []
    for archive_name, our_name, label in STATIONS:
        ref_path = os.path.join(ARCHIVE_DIR, archive_name)
        our_path = os.path.join(args.run_dir, our_name)
        if not os.path.isfile(ref_path):
            print(f'  {label}: MISSING reference file {ref_path} -- skipped')
            continue
        if not os.path.isfile(our_path):
            print(f'  {label}: MISSING run file {our_path} -- skipped')
            continue
        r = compare_one(ref_path, our_path, label)
        results.append(r)
        print(f'  {label}')
        print(f'    ref  : {r["n_rows_ref"]} rows, t in {r["t_span_ref"]}')
        print(f'    ours : {r["n_rows_ours"]} rows, t in {r["t_span_ours"]}')
        for col, v in r['maxdiff_by_column'].items():
            print(f'    max|diff| {col:14s}: {v:.4e}')

    print('\nThis script is REPORT-ONLY (module docstring): it never asserts '
          'and always exits 0. No scalar bound is recorded yet -- re-run '
          'once a spec-resolution (100 m, 8 s) test.tpv12 run exists.')

    if not args.no_json:
        os.makedirs(OUT_DIR, exist_ok=True)
        ts = datetime.datetime.utcnow().strftime('%Y%m%dT%H%M%SZ')
        out_path = os.path.join(OUT_DIR, f'evidence_tpv12_{ts}.json')
        with open(out_path, 'w') as f:
            json.dump(dict(provenance=prov, case='tpv12', run_dir=args.run_dir,
                            gate_resolution=args.gate_resolution,
                            rupture_counts=rc, results=results),
                      f, indent=2, default=str)
        print(f'\nWrote {out_path}')
    return 0


if __name__ == '__main__':
    sys.exit(main())

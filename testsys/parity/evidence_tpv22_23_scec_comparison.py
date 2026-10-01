#! /usr/bin/env python3
"""
Independent cross-code validation for TPV22/TPV23 (rule 17 step 6), the same
kind of check test.tpv29/test.tpv30 have (testsys/parity/
evidence_tpv29_scec_comparison.py) -- a committed SCRIPT, not prose-only with
a gitignored baseline (TPV29's original 2015 write-up did that; pathway item
28; rule 17 step 6 names it as the thing to NOT repeat).

REPORT-ONLY, same contract as evidence_tpv29_scec_comparison.py: never a
gate, never wired into testsys/run.py, invoked by hand. Prints the comparison
and says REPRODUCES/DRIFTED/EXPECTED-MISMATCH per metric; it never asserts
and always exits 0.

WHAT IS COMPARED AND WHERE IT COMES FROM
-----------------------------------------
Two INDEPENDENT codes' public SCEC cvws submissions for TPV22/TPV23, fetched
2026-10-01 from the portal's "Public Area" (https://strike.scec.org/cvws/
cgi-bin/cvws.cgi, the same CGI scripts/scec/fetch_all_dliu.py automates for
this repo's OWN submissions -- this script's reference data is OTHER
submitters', so it was fetched by hand, see scratch/specs/tpv2223/
scec_cross_code/<case>/<user>/PROVENANCE, if present, or this module's own
citation below) against a completed run of EQdyna's current code:

  1. THE REFERENCE -- scratch/specs/tpv2223/scec_cross_code/<case>/<user>/,
     <case> in (tpv22, tpv23), <user> in:
       kaneko -- Yoshihiro Kaneko, SPECFEM3D (spectral element), 100 m,
                 fetched 2026-10-01 (cvws benchmark='tpv22'/'tpv23', user
                 'kaneko', files cplot_1/cplot_2/fault1st000dp000/
                 fault2st000dp000/fault2st050dp050). INDEPENDENT of EQdyna
                 (different method, different group) -- the primary
                 cross-code check.
       payne  -- Ryan Payne (submitter)/Benchun Duan (author, per the
                 file's own header), EQdyna (finite element), 100 m, same
                 fetch. This is EQdyna ITSELF (an earlier run, 2013), so it
                 is NOT an independent cross-code check -- kept alongside
                 kaneko only as a same-code sanity cross-check (does
                 EQdyna's own 2013 submission agree with EQdyna's own
                 2015-vintage SPECFEM3D peer at the stations both ship?),
                 reported separately and never averaged with kaneko's.
     Format: 8-column on-fault time series (t, h-slip, h-slip-rate,
     h-shear-stress, v-slip, v-slip-rate, v-shear-stress, n-stress) -- the
     SCEC standard TPV22/23 on-fault format (TPV22_23_Description_v08, Part
     6, p.15), identical column layout to this repo's own faultst*.txt, so
     no reformatting is needed, only a time-window intersection (see below).
     fault1st000dp000 = fault #1, 0 km along-strike, 0 km down-dip.
     fault2st000dp000 = fault #2, 0 km along-strike, 0 km down-dip (the
     node geometrically closest to fault #1's overlap zone).
     fault2st050dp050 = fault #2, 5 km along-strike, 5 km down-dip (the
     station that showed the extensional/compressional stress-sign split
     in this mission's own gate sanity check -- see case_input/test.tpv22's
     README/the mission report for the measured numbers).

  2. THE CURRENT-CODE RUN -- a completed test.tpv22/test.tpv23 case
     directory (NOT launched here -- rule 9, same as evidence_tpv29's
     policy: case setup + a real solver run is not "cheap"). Pass
     --run-dir /path/to/completed/case. This script does NOT require the
     fortran binary or a fresh run; it only reads files already on disk.

WHY NO BOUND IS ASSERTED YET (read before expecting PASS/FAIL)
----------------------------------------------------------------
This repo's own tpv22/tpv23 gate (case_input/test.tpv22, test.tpv23) runs at
GATE-COARSE resolution (dx=dz=1000 m, dy=400/500 m, 5 s term -- see
tpv22_23_common.py) -- a REGRESSION check, not a claim about physical
accuracy (same framing scripts/scec/README.md gives tpv29/tpv30's coarse
gate cells). The archive's kaneko/payne submissions are at the SPEC
resolution (100 m) and the full 15 s duration (full_specs.py's test.tpv22/
test.tpv23 entries, NOT run by this mission -- rule 17 step 5, scheduling
only). Comparing a 1000 m/5 s run against a 100 m/15 s reference will show
LARGE differences almost everywhere past the first second or two -- that is
an EXPECTED RESOLUTION MISMATCH, not evidence of a bug, and this script
reports it as such rather than asserting a bound no coarse run could pass.
Once a spec-resolution (100 m, 15 s) run exists, re-running this script
against that run directory is what makes the comparison meaningful as a
pass/fail check; a RECORDED number to check future runs against belongs
here only after that run exists (rule 4 -- no claim before a measurement).

Run:
    python3 testsys/parity/evidence_tpv22_23_scec_comparison.py --case tpv22 --run-dir /path/to/completed/test.tpv22/case
    python3 testsys/parity/evidence_tpv22_23_scec_comparison.py --case tpv23 --run-dir /path/to/completed/test.tpv23/case --no-json

Writes a timestamped JSON snapshot to testsys/parity/evidence_output/
(gitignored -- a number tied to whatever run directory was passed is
evidence-of-a-run, not golden reference data; rule 4/rule 7).
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
ARCHIVE_ROOT = os.path.join(REPO_ROOT, 'scratch', 'specs', 'tpv2223', 'scec_cross_code')

# Our own case's file name for each archive station, per case (fault #1 has
# no 'ft'-prefix at ntotft>1's fault-1 convention; fault #2 uses the
# faultTag() prefix 'ft2_' -- scripts/lib.py:faultTag, case_input/test.tpv22/
# tpv22_23_common.py's st_coor_on_fault_per_fault).
STATIONS = [
    # (archive file name, our file name, human label)
    ('fault1st000dp000.txt', 'faultst000dp000.txt',
     'fault #1, 0 km strike, 0 km down-dip'),
    ('fault2st000dp000.txt', 'faultstft2_000dp000.txt',
     'fault #2, 0 km strike, 0 km down-dip (overlap-zone edge)'),
    ('fault2st050dp050.txt', 'faultstft2_050dp050.txt',
     'fault #2, 5 km strike, 5 km down-dip'),
]
SOURCES = [
    ('kaneko', 'Yoshihiro Kaneko, SPECFEM3D (spectral element), 100 m, '
               '2013 -- INDEPENDENT of EQdyna'),
    ('payne', 'Ryan Payne / Benchun Duan, EQdyna (finite element), 100 m, '
              '2013 -- SAME CODE FAMILY, sanity cross-check only, NOT '
              'independent'),
]


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
    """For every row of `ours` (assumed the SHORTER/coarser series), find
    the nearest-time row of `ref` and return (ours_rows, matched_ref_rows)
    aligned 1:1 -- a resampling-free comparison (never interpolates a
    value), same spirit as evidence_tpv29's exact-grid-lookup approach."""
    t_ref = ref[:, 0]
    out_ref = np.empty((ours.shape[0], ref.shape[1]))
    for i, t in enumerate(ours[:, 0]):
        j = int(np.argmin(np.abs(t_ref - t)))
        out_ref[i] = ref[j]
    return out_ref


COLS = ['h-slip (m)', 'h-slip-rate (m/s)', 'h-shear-stress (MPa)',
        'v-slip (m)', 'v-slip-rate (m/s)', 'v-shear-stress (MPa)',
        'n-stress (MPa)']


def compare_one(ref_path, our_path, label):
    ref = _numeric_rows(ref_path)
    ours = _numeric_rows(our_path)
    if ref.shape[1] != 8 or ours.shape[1] != 8:
        raise SystemExit(f'{label}: expected 8 columns, got ref={ref.shape[1]} '
                         f'ours={ours.shape[1]}')
    t_span_ref = (float(ref[0, 0]), float(ref[-1, 0]))
    t_span_ours = (float(ours[0, 0]), float(ours[-1, 0]))
    matched_ref = _nearest_time_match(ref, ours)
    diff = np.abs(ours[:, 1:] - matched_ref[:, 1:])
    maxdiff = diff.max(axis=0)
    # Rupture onset: first time |h-slip| exceeds 1 mm, in each series' own
    # time base (not matched -- this is each series' own arrival time).
    def onset(a):
        idx = np.nonzero(np.abs(a[:, 1]) > 1.0e-3)[0]
        return float(a[idx[0], 0]) if idx.size else None
    return dict(
        label=label, t_span_ref=t_span_ref, t_span_ours=t_span_ours,
        n_rows_ref=int(ref.shape[0]), n_rows_ours=int(ours.shape[0]),
        maxdiff_by_column=dict(zip(COLS, (float(x) for x in maxdiff))),
        onset_ref_s=onset(ref), onset_ours_s=onset(ours),
    )


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
    ap.add_argument('--case', required=True, choices=('tpv22', 'tpv23'))
    ap.add_argument('--run-dir', required=True,
                     help='completed test.tpv22/test.tpv23 case directory '
                          '(any resolution/term -- a resolution mismatch '
                          'against the 100 m archive is reported, not hidden)')
    ap.add_argument('--no-json', action='store_true')
    args = ap.parse_args()

    archive_dir = os.path.join(ARCHIVE_ROOT, args.case)
    if not os.path.isdir(archive_dir):
        raise SystemExit(f'{archive_dir}: missing -- fetch per this module\'s '
                         f'docstring (scratch/specs/tpv2223/scec_cross_code/)')

    prov = provenance()
    print('==== provenance ====')
    for k, v in prov.items():
        print(f'  {k:9s}: {v}')
    print(f'\ncase: {args.case}')
    print(f'run_dir: {args.run_dir}')

    results = {}
    for source, source_label in SOURCES:
        print(f'\n==== source: {source} ({source_label}) ====')
        results[source] = []
        for archive_name, our_name, label in STATIONS:
            ref_path = os.path.join(archive_dir, source, archive_name)
            our_path = os.path.join(args.run_dir, our_name)
            if not os.path.isfile(ref_path):
                print(f'  {label}: MISSING reference file {ref_path} -- skipped')
                continue
            if not os.path.isfile(our_path):
                print(f'  {label}: MISSING run file {our_path} (station not '
                      f'produced at this run\'s resolution?) -- skipped')
                continue
            r = compare_one(ref_path, our_path, label)
            results[source].append(r)
            same_window = abs(r['t_span_ours'][1] - r['t_span_ref'][1]) < 1.0
            tag = 'COMPARABLE (similar time windows)' if same_window else \
                  'EXPECTED-MISMATCH (resolution/term differ -- see module docstring)'
            print(f'  {label}')
            print(f'    ref  : {r["n_rows_ref"]} rows, t in {r["t_span_ref"]}')
            print(f'    ours : {r["n_rows_ours"]} rows, t in {r["t_span_ours"]} -- {tag}')
            print(f'    rupture onset (|h-slip|>1mm): ref={r["onset_ref_s"]}, ours={r["onset_ours_s"]}')
            for col, v in r['maxdiff_by_column'].items():
                print(f'    max|diff| {col:22s}: {v:.4e}')

    print('\nThis script is REPORT-ONLY (module docstring): it never asserts '
          'and always exits 0. No scalar bound is recorded yet -- see "WHY NO '
          'BOUND IS ASSERTED YET" above. Re-run once a spec-resolution '
          '(100 m, 15 s) test.tpv22/test.tpv23 run exists to get a '
          'comparison that is not dominated by the resolution/term mismatch.')

    if not args.no_json:
        os.makedirs(OUT_DIR, exist_ok=True)
        ts = datetime.datetime.utcnow().strftime('%Y%m%dT%H%M%SZ')
        out_path = os.path.join(OUT_DIR, f'evidence_{args.case}_{ts}.json')
        with open(out_path, 'w') as f:
            json.dump(dict(provenance=prov, case=args.case, run_dir=args.run_dir,
                            results=results), f, indent=2, default=str)
        print(f'\nWrote {out_path}')
    return 0


if __name__ == '__main__':
    sys.exit(main())

#! /usr/bin/env python3
"""
Compares EQdyna's OWN committed TPV30 Fortran reference
(test.reference.results/test.tpv30/, dx=500 m, rule 17 step 4 gate run) against
the OWNER'S OWN published TPV30 submission on the public SCEC/USGS cvws
archive (scec_archive/tpv30/eqdyna-v3.1-100m-2015/, dx=100 m) -- rule 17 step
6's independent-validation requirement, previously unsatisfied for TPV30
(NOTES_tpv30_gate.md).

REPORT-ONLY, same pattern as evidence_tpv29_scec_comparison.py and
evidence_tpv30_vs_tpv29_contrast.py in this directory: never a gate, never
wired into testsys/run.py, invoked by hand. Never asserts, always exits 0.

WHAT THIS SCRIPT DOES AND DOES NOT CLAIM
-----------------------------------------
This is a CROSS-RESOLUTION comparison (500 m current-code gate run vs 100 m
2015 archive submission), NOT a same-resolution parity check like
evidence_tpv29_scec_comparison.py's 100-vs-100 comparison. A 500 m run is
expected to disagree with a 100 m run by a lot -- the owner's own 100 m
submission itself ranks 14/14 (farthest from the group median) among 14
independent TPV30 submissions at the SCEC portal, and rank improves
monotonically 100m -> 50m -> 25m in the owner's own record (100m:14/14,
50m:12/14, 25m:8/14, per scec_archive/tpv30/*/PROVENANCE.md). So this script
answers ONE question only: is EQdyna's 500 m run RECOGNIZABLY THE SAME
PHYSICS as the owner's own finer run (same hypocenter, same rupture
direction/extent, same order-of-magnitude slip and moment), not "does it
match to any accuracy standard". Treat every DRIFTED verdict below as a
regression flag, and every REPRODUCES verdict as "recognizably the same
physics", never as an accuracy claim. See NOTES_tpv30_gate.md Step 2 for the
mission framing this script exists to answer.

This script does NOT touch the separate, open, unresolved finding that
python-numpy/python-jax diverge from EQdyna's OWN Fortran by up to 30% at
t=20s on this same case (NOTES_tpv30_gate.md, "STOP -- finding, not a
landing"). That is a PORT-CORRECTNESS question (Fortran vs Python) already
diagnosed by binary search in that document. This script only ever compares
Fortran-vs-Fortran (EQdyna vs EQdyna, 21 years apart) -- a PHYSICS-VALIDITY
question. The two questions are independent; see this script's own
print_report for both verdicts stated side by side, never conflated.

WHERE EACH INPUT COMES FROM
----------------------------
  1. THE 2015 100 m SUBMISSION -- scec_archive/tpv30/eqdyna-v3.1-100m-2015/
     (same CGI fetch path/format as TPV29's; see that directory's own
     PROVENANCE.md). `cplot`: 3-column strike/dip/rupture-time grid, 401x201
     nodes, sentinel 1000.0 s -- confirmed here to be the SAME 40kmx20km rough
     fault frame as TPV29's own 100 m archive (same node count, same extent),
     consistent with the spec's "TPV30 = TPV29 + material change only".
     24 faultst<S>dp<D> on-fault station time series ship alongside it.
  2. EQdyna's OWN 500 m REFERENCE --
     test.reference.results/test.tpv30/frt.canonical.txt, the ALREADY
     deduped+lexsorted 22-column fault-node array frozen by the rule-17-step-4
     gate run (NOTES_tpv30_gate.md). This is a SINGLE static file, not a live
     run directory -- there is no separate SCECRuptureTime.txt/faultst*.txt
     for this coarse gate run, so this script reads final rupture time and
     final slip straight out of the canonical array's own columns (x=col0,
     z=col2, rupture_time=col3, slipStrike=col4, slipDip=col5 -- confirmed
     against testsys/frt_canonical.py's own column layout and
     evidence_tpv30_vs_tpv29_contrast.py's identical column use).
     dx/fault-extent/hypocenter are read from case_input/test.tpv30/
     user_defined_params.py directly (NOT from testsys/matrix.py, which no
     longer carries a test.tpv30 entry -- it was deliberately removed when
     the case was found not ready to gate; see NOTES_tpv30_gate.md).

COORDINATE MAPPING (verified, not assumed)
-------------------------------------------
frt.canonical.txt's z column runs -20000..0 (0 = free surface, matching
par.fzmin/fzmax = -20000/0); the archive's cplot "k" (dip) column runs
0..20000 (0 = free surface). So dip_m = -z. x (strike) is directly comparable,
no sign flip -- both run -20000..20000. Checked directly against the
reference file's own x/z ranges before trusting this mapping (main() smoke
check, raises rather than silently proceeding if it disagrees).

Because dx=500 exactly divides the archive's dx=100, EVERY ONE of the current
run's 3321 fault nodes has an EXACT coordinate match in the 100 m archive grid
-- this script does a direct dict lookup, never nearest-neighbour or
interpolation, for the full-grid rupture-time comparison (METRIC 1). It
RAISES if any node fails to find an exact match (a grid-alignment bug would
show up as a hole in the comparison, not as noise). Slip comparison (METRIC 3)
is different: the archive ships only 24 named station time series (no full 2D
slip field, same limitation evidence_tpv29_scec_comparison.py documents for
TPV29), and most of those station coordinates are NOT on the 500 m current
grid (e.g. faultst-042dp061 = -4200 m, not a multiple of 500). For those,
this script falls back to the current run's NEAREST 500 m grid node, and
reports the offset distance next to each number -- always labelled
"nearest-node", never claimed to be the same physical point exactly.

Run:
    python3 testsys/parity/evidence_tpv30_scec_comparison.py
    python3 testsys/parity/evidence_tpv30_scec_comparison.py --submission-res 50m
    python3 testsys/parity/evidence_tpv30_scec_comparison.py --no-json

Writes a timestamped JSON snapshot to testsys/parity/evidence_output/
(gitignored -- rule 4/rule 7: a number tied to a specific archive resolution
choice is evidence-of-a-comparison, not golden reference data).
"""
import argparse
import datetime
import glob
import os
import socket
import subprocess
import sys
import json

import numpy as np

TESTSYS = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(TESTSYS))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)
OUT_DIR = os.path.join(TESTSYS, 'evidence_output')

# Reuse the generic archive-file I/O this script shares with TPV29's own
# comparison (rule 1 -- do not duplicate): the archive file FORMATS are
# byte-identical between the two benchmarks (same CGI, same author, same
# fault frame), only the DATA differs.
from testsys.parity import evidence_tpv29_scec_comparison as tpv29cmp  # noqa: E402

CASE = 'test.tpv30'
REFERENCE_FILE = os.path.join(
    REPO_ROOT, 'test.reference.results', CASE, 'frt.canonical.txt')
# 100m is the ONLY resolution this script's --submission-res actually
# supports: tpv29cmp.load_rupture_grid hardcodes NX_EXPECT=401/NZ_EXPECT=201/
# DX_M=100 as MODULE constants (shared, correctly, with TPV29's own 100 m
# comparison), so it raises loudly on the 50 m (801x401) and 25 m (1601x801)
# archives rather than silently misreading them -- confirmed by actually
# running --submission-res 50m/25m against real files, not assumed. Making
# those two work would mean generalizing that shared function's grid-shape
# check, which is out of scope here (the mission's own Step 4 asks about
# COSTING OUT a fresh 100 m EQdyna run, not about reading the 50/25 m
# archives) and not done. The 50 m/25 m paths are recorded for provenance
# only; SUPPORTED_SUBMISSION_RES below is what --choices actually allows.
SUBMISSION_DIRS = {
    '100m': os.path.join(REPO_ROOT, 'scec_archive', 'tpv30', 'eqdyna-v3.1-100m-2015'),
    '50m': os.path.join(REPO_ROOT, 'scec_archive', 'tpv30', 'eqdyna-v3.1-50m-2015'),
    '25m': os.path.join(REPO_ROOT, 'scec_archive', 'tpv30', 'eqdyna-v3.1-25m-2015'),
}
SUPPORTED_SUBMISSION_RES = ('100m',)

RHO, VS = tpv29cmp.RHO, tpv29cmp.VS   # TPV30 spec p.11: same material as TPV29
SENTINEL_CURRENT = 9.0e4              # EQdyna's own "never ruptured" sentinel
SENTINEL_ARCHIVE = tpv29cmp.UNRUPTURED  # 900, separates < 1000.0 s archive sentinel


def _case_params(case):
    """par.dx, fault extent, and hypocenter read straight from the case's own
    compset -- mirrors evidence_tpv30_vs_tpv29_contrast.py's _dx_and_nproc /
    hypocenter helpers, but WITHOUT the matrix.FORTRAN_RANKS lookup those use
    (test.tpv30 has no matrix.py entry -- it was deliberately removed,
    NOTES_tpv30_gate.md -- and this script never needs a rank count, since it
    reads the single already-canonicalized frt.canonical.txt, not per-rank
    frt.txt<rank> files)."""
    case_input = os.path.join(REPO_ROOT, 'case_input', case)
    scripts_dir = os.path.join(REPO_ROOT, 'scripts')
    # user_defined_params.py does `from defaultParameters import *` /
    # `from lib import *` with no sys.path handling of its own -- it is
    # normally imported by scripts/case.setup, which lives IN scripts/ and so
    # gets that directory on sys.path for free (Python's own script-dir
    # rule). Standalone here, scripts/ has to be added explicitly the same
    # way case_input does -- this is a pre-existing gap shared by
    # evidence_tpv30_vs_tpv29_contrast.py's identical _dx_and_nproc/
    # hypocenter helpers (untested to completion there per
    # NOTES_tpv30_gate.md; caught here by actually running this script).
    sys.path.insert(0, scripts_dir)
    sys.path.insert(0, case_input)
    cwd = os.getcwd()
    try:
        saved = {k: sys.modules.pop(k) for k in
                 ('user_defined_params', 'defaultParameters', 'lib')
                 if k in sys.modules}
        # tpv29GeometryTools.loadEQdynaGeometry() opens
        # 'bFault_Rough_Geometry.tpv29.100m.txt' by a bare relative name
        # (same file case.setup's own cwd convention expects: it always runs
        # from inside the case directory) -- chdir there for the duration of
        # the import, same reason scripts_dir is added above.
        os.chdir(case_input)
        from user_defined_params import par  # noqa
        out = dict(dx=float(par.dx),
                   fxmin=float(par.fxmin), fxmax=float(par.fxmax),
                   fzmin=float(par.fzmin), fzmax=float(par.fzmax),
                   xsource=float(par.xsource), zsource=float(par.zsource))
    finally:
        os.chdir(cwd)
        sys.path.remove(case_input)
        sys.path.remove(scripts_dir)
        for k, v in saved.items():
            sys.modules[k] = v
        sys.modules.pop('user_defined_params', None)
    return out


def load_reference(path, p):
    """frt.canonical.txt -> dict of (nz,nx) fields, same layout as
    evidence_tpv30_vs_tpv29_contrast.load_trial, but from ONE
    already-canonicalized file instead of a live run directory's
    frt.txt<rank> glob (no dedupe needed -- frt_canonical.py already did it,
    see testsys/frt_canonical.py's own docstring)."""
    if not os.path.isfile(path):
        raise SystemExit(
            f'{path}: no such file. This is the committed rule-17-step-4 gate '
            f'reference (NOTES_tpv30_gate.md) -- it should already be in the '
            f'repo; do not regenerate it from a fresh run for this comparison '
            f'(that would be a different, undocumented run, rule 4).')
    a = np.loadtxt(path)
    if a.shape[1] != 22:
        raise SystemExit(f'{path}: {a.shape[1]} columns, expected 22 '
                          f'(testsys/frt_canonical.py FRT_COLUMNS)')
    dx = p['dx']
    nx = round((p['fxmax'] - p['fxmin']) / dx) + 1
    nz = round((p['fzmax'] - p['fzmin']) / dx) + 1
    x = p['fxmin'] + np.arange(nx) * dx
    z = p['fzmin'] + np.arange(nz) * dx
    rupt = np.full((nz, nx), np.nan)
    slipS = np.full((nz, nx), np.nan)
    slipD = np.full((nz, nx), np.nan)
    ii = np.rint((a[:, 0] - p['fxmin']) / dx).astype(int)
    jj = np.rint((a[:, 2] - p['fzmin']) / dx).astype(int)
    if a.shape[0] != nx * nz:
        raise SystemExit(f'{path}: {a.shape[0]} rows, expected {nx * nz} '
                          f'({nx}x{nz} at dx={dx:g} over the case\'s own '
                          f'fault extent) -- wrong file or wrong case params')
    rupt[jj, ii] = a[:, 3]
    slipS[jj, ii] = a[:, 4]
    slipD[jj, ii] = a[:, 5]
    if np.isnan(rupt).any():
        raise SystemExit(f'{path}: {int(np.isnan(rupt).sum())} nodes missing '
                          f'after gridding -- coordinate mapping is wrong, '
                          f'not just coarse-vs-fine noise')
    return dict(x=x, z=z, dip=-z, rupt=rupt, slipS=slipS, slipD=slipD, dx=dx)


def metric1_rupture_time_full_grid(ref, x15, z15, T15):
    """Every current-run node has an EXACT coordinate match in the 100 m
    archive grid (dx=500 divides dx=100 exactly) -- direct dict lookup, never
    nearest-neighbour. Raises if a node is missing (grid misalignment bug),
    per this module's own coordinate-mapping section."""
    lut15 = {(round(x15[j], 1), round(z15[i], 1)): T15[i, j]
             for i in range(z15.size) for j in range(x15.size)}
    nz, nx = ref['rupt'].shape
    dt, missing = [], []
    both_ruptured = both_never = only_current = only_archive = 0
    for iz in range(nz):
        for ix in range(nx):
            key = (round(float(ref['x'][ix]), 1), round(float(ref['dip'][iz]), 1))
            if key not in lut15:
                missing.append(key)
                continue
            t15 = lut15[key]
            tnow = ref['rupt'][iz, ix]
            r15 = t15 < SENTINEL_ARCHIVE
            rnow = tnow < SENTINEL_CURRENT
            if r15 and rnow:
                both_ruptured += 1
                dt.append(abs(t15 - tnow))
            elif (not r15) and (not rnow):
                both_never += 1
            elif rnow and not r15:
                only_current += 1
            else:
                only_archive += 1
    if missing:
        raise SystemExit(
            f'{len(missing)} of {nx * nz} current-run nodes have NO exact '
            f'coordinate match in the 100 m archive grid (first: {missing[0]}) '
            f'-- coordinate mapping is wrong, this is not expected for a '
            f'dx=500-divides-dx=100 grid.')
    return dict(n_nodes=nx * nz, n_both_ruptured=both_ruptured,
                n_both_never=both_never, n_only_current=only_current,
                n_only_archive=only_archive,
                median_dt_s=float(np.median(dt)) if dt else None,
                mean_dt_s=float(np.mean(dt)) if dt else None,
                max_dt_s=float(np.max(dt)) if dt else None)


def metric2_ruptured_area(ref, x15, z15, T15):
    """Corner + naive area rules (evidence_tpv29_scec_comparison's own
    formulas, reused not reimplemented), computed independently on each
    grid's own resolution -- NOT expected to agree closely: a coarser mesh's
    nucleation patch and rupture front are quantized to a 500 m cell, the
    archive's to a 100 m cell, so discretization noise alone gives a
    real difference even under identical physics.

    rupture_area_km2/naive_ruptured_area_km2 take dz = np.diff(z).mean(),
    which is SIGNED -- ref['dip'] runs 20000 -> 0 (descending, because it is
    built as -ref['z'] and ref['z'] itself runs -20000 -> 0 ascending), so it
    must be flipped to ascending order here before calling, or the area comes
    back negative (caught by actually running this script against real data,
    not assumed)."""
    dip_asc = ref['dip'][::-1]
    rupt_asc = ref['rupt'][::-1, :]
    area_now = tpv29cmp.rupture_area_km2(
        ref['x'], dip_asc, rupt_asc, sentinel=SENTINEL_CURRENT)
    naive_now = tpv29cmp.naive_ruptured_area_km2(
        ref['x'], dip_asc, rupt_asc, sentinel=SENTINEL_CURRENT)
    area_15 = tpv29cmp.rupture_area_km2(x15, z15, T15, sentinel=SENTINEL_ARCHIVE)
    naive_15 = tpv29cmp.naive_ruptured_area_km2(x15, z15, T15, sentinel=SENTINEL_ARCHIVE)
    return dict(area_now_km2=area_now, area_2015_km2=area_15,
                naive_now_km2=naive_now, naive_2015_km2=naive_15)


def metric3_final_slip_nearest_node(ref, submission_dir):
    """24 archive stations vs the current run's NEAREST 500 m grid node
    (labelled, never claimed exact -- see module docstring)."""
    stations = sorted(glob.glob(os.path.join(submission_dir, 'faultst*')))
    rows = []
    for path in stations:
        name = os.path.basename(path)
        strike_m, dip_m = tpv29cmp.station_coords_m(name)
        ix = int(np.argmin(np.abs(ref['x'] - strike_m)))
        iz = int(np.argmin(np.abs(ref['dip'] - dip_m)))
        offset_m = float(np.hypot(ref['x'][ix] - strike_m, ref['dip'][iz] - dip_m))
        slip_archive = tpv29cmp.station_final_slip(path)
        slip_now = float(np.hypot(ref['slipS'][iz, ix], ref['slipD'][iz, ix]))
        pct = abs(slip_now - slip_archive) / slip_archive * 100.0 if slip_archive else None
        rows.append(dict(station=name, strike_m=strike_m, dip_m=dip_m,
                          nearest_node_offset_m=offset_m,
                          slip_2015_m=slip_archive, slip_now_m=slip_now,
                          slip_pct_diff=pct))
    pcts = [r['slip_pct_diff'] for r in rows if r['slip_pct_diff'] is not None]
    return dict(per_station=rows,
                n_stations=len(rows),
                median_slip_pct=float(np.median(pcts)) if pcts else None,
                max_nearest_offset_m=float(max(r['nearest_node_offset_m'] for r in rows))
                if rows else None)


def metric4_moment_mw_current(ref):
    """Current (500 m) run only -- the archive has no full 2D slip field, same
    limitation evidence_tpv29_scec_comparison.load_frt_moment documents.
    Same formula (M0 = mu*sum(slip*dx*dz), Hanks-Kanamori Mw), reimplemented
    directly over the already-gridded (nz,nx) field here rather than through
    that function, because this run's input is one pre-canonicalized file, not
    a frt.txt<rank> glob to dedupe."""
    mu = RHO * VS ** 2
    slip = np.hypot(ref['slipS'], ref['slipD'])
    moment = float(np.nansum(slip) * ref['dx'] * ref['dx'] * mu)
    mw = 2.0 / 3.0 * np.log10(moment * 1.0e7) - 10.7
    return dict(moment_Nm=moment, mw=float(mw), n_nodes=int(slip.size))


def provenance():
    sha = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=REPO_ROOT,
                          capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(['git', 'status', '--porcelain'], cwd=REPO_ROOT,
                                 capture_output=True, text=True).stdout.strip())
    return dict(date_utc=datetime.datetime.utcnow().isoformat() + 'Z',
                host=socket.gethostname(), sha=sha, dirty=dirty)


def print_report(res, submission_res, submission_dir):
    m1, m2, m3, m4 = res['metric1'], res['metric2'], res['metric3'], res['metric4']
    print(f'\n=== TPV30: EQdyna 500 m gate reference vs EQdyna\'s own '
          f'{submission_res} 2015 SCEC submission ===')
    print(f'reference : {REFERENCE_FILE}')
    print(f'submission: {submission_dir}')
    print('\nFraming (see module docstring): this is a REGRESSION check '
          '("recognizably the same physics"), never an accuracy claim -- a '
          '500 m run is expected to disagree substantially with a 100 m run.')

    print('\n==== METRIC 1: rupture-time diff, FULL 500 m grid, exact-coordinate '
          'match against the 100 m archive (%d nodes) ====' % m1['n_nodes'])
    print(f'  both ruptured           : {m1["n_both_ruptured"]}')
    print(f'  both never ruptured     : {m1["n_both_never"]}')
    print(f'  ruptured now only       : {m1["n_only_current"]}')
    print(f'  ruptured in archive only: {m1["n_only_archive"]}')
    if m1['median_dt_s'] is not None:
        print(f'  |t_500m - t_2015_100m|: median {m1["median_dt_s"]:.3f} s, '
              f'mean {m1["mean_dt_s"]:.3f} s, max {m1["max_dt_s"]:.3f} s')
    extent_ratio = m1['n_both_ruptured'] / max(1, m1['n_both_ruptured'] + m1['n_only_archive'])
    print(f'  -> {extent_ratio*100:.1f}% of archive-ruptured nodes also ruptured '
          f'at 500 m (extent agreement, not timing accuracy)')

    print('\n==== METRIC 2: ruptured area (each grid at its own resolution) ====')
    print(f'  corner rule : 500m now = {m2["area_now_km2"]:.1f} km2   '
          f'2015 100m = {m2["area_2015_km2"]:.1f} km2   '
          f'ratio = {m2["area_now_km2"]/m2["area_2015_km2"]:.3f}')
    print(f'  naive rule  : 500m now = {m2["naive_now_km2"]:.1f} km2   '
          f'2015 100m = {m2["naive_2015_km2"]:.1f} km2   '
          f'ratio = {m2["naive_now_km2"]/m2["naive_2015_km2"]:.3f}')

    print(f'\n==== METRIC 3: final slip, {m3["n_stations"]} archive stations vs '
          f'NEAREST 500 m current-run node (max offset '
          f'{m3["max_nearest_offset_m"]:.0f} m) ====')
    print(f'  median |slip_now - slip_2015| / slip_2015 = '
          f'{m3["median_slip_pct"]:.1f}%')
    for r in m3['per_station']:
        print(f'    {r["station"]:20s} 2015={r["slip_2015_m"]:6.3f} m  '
              f'now(nearest,{r["nearest_node_offset_m"]:5.0f} m away)='
              f'{r["slip_now_m"]:6.3f} m  diff={r["slip_pct_diff"]:6.1f}%')

    print('\n==== METRIC 4: Mw (500 m current run only -- archive has no full '
          '2D slip field, see docstring) ====')
    print(f'  M0 = {m4["moment_Nm"]:.4e} N*m over {m4["n_nodes"]} nodes, '
          f'rho={RHO:g} kg/m3, Vs={VS:g} m/s')
    print(f'  Mw = {m4["mw"]:.3f}')

    print('\n==== Verdict (Step 2/3 framing, NOTES_tpv30_gate.md) ====')
    same_physics = (extent_ratio > 0.85 and
                     0.3 < m2['area_now_km2'] / m2['area_2015_km2'] < 3.0 and
                     m4['moment_Nm'] > 0)
    print(f'  (A) physics validity (this script): '
          f'{"ROUGHLY MATCHES" if same_physics else "DOES NOT ROUGHLY MATCH"} '
          f'the owner\'s own {submission_res} published TPV30 submission -- '
          f'{extent_ratio*100:.0f}% rupture-extent overlap, '
          f'area ratio {m2["area_now_km2"]/m2["area_2015_km2"]:.2f}, '
          f'Mw {m4["mw"]:.2f}.')
    print('  (B) port correctness (numpy/jax vs Fortran): UNCHANGED by this '
          'script -- NOTES_tpv30_gate.md\'s prior finding stands (Fortran vs '
          'python-numpy/python-jax diverge up to 30% at t=20s, root cause not '
          'found). This script never touches that question.')
    print('  These are two SEPARATE questions (mission Step 3) -- (A) matching '
          'does NOT imply (B) is resolved, and no gate action follows from '
          'this script by itself.')

    print('\nThis script is REPORT-ONLY (module docstring): it never asserts '
          'and always exits 0.')
    return same_physics


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--submission-res', default='100m', choices=SUPPORTED_SUBMISSION_RES,
                     help='which archived resolution to compare against '
                          '(only 100m is supported -- see SUBMISSION_DIRS\' '
                          'own comment for why 50m/25m are not; 100m is also '
                          'the mission\'s own choice, the cheapest resolution '
                          'giving a real cross-code point)')
    ap.add_argument('--no-json', action='store_true')
    args = ap.parse_args()

    submission_dir = SUBMISSION_DIRS[args.submission_res]
    if not os.path.isfile(os.path.join(submission_dir, 'cplot')):
        raise SystemExit(f'{submission_dir}: no cplot file')

    p = _case_params(CASE)
    ref = load_reference(REFERENCE_FILE, p)
    x15, z15, T15 = tpv29cmp.load_rupture_grid(os.path.join(submission_dir, 'cplot'))

    # Smoke-check the coordinate mapping claimed in the module docstring
    # BEFORE trusting any metric built on it.
    if not (float(x15.min()) == float(ref['x'].min()) == p['fxmin']
            and float(x15.max()) == float(ref['x'].max()) == p['fxmax']):
        raise SystemExit('strike (x) extents disagree between the 500 m '
                          'reference and the archive grid -- coordinate '
                          'mapping assumption is wrong, stopping before '
                          'computing misleading metrics.')
    if not (float(z15.min()) == 0.0 and float(z15.max()) == -p['fzmin']
            and float(ref['dip'].min()) == 0.0 and float(ref['dip'].max()) == -p['fzmin']):
        raise SystemExit('dip extents disagree between the 500 m reference '
                          '(dip = -z) and the archive grid -- coordinate '
                          'mapping assumption is wrong, stopping before '
                          'computing misleading metrics.')

    prov = provenance()
    print('==== provenance ====')
    for k, v in prov.items():
        print(f'  {k:9s}: {v}')

    res = dict(
        metric1=metric1_rupture_time_full_grid(ref, x15, z15, T15),
        metric2=metric2_ruptured_area(ref, x15, z15, T15),
        metric3=metric3_final_slip_nearest_node(ref, submission_dir),
        metric4=metric4_moment_mw_current(ref),
    )
    same_physics = print_report(res, args.submission_res, submission_dir)

    if not args.no_json:
        os.makedirs(OUT_DIR, exist_ok=True)
        ts = datetime.datetime.utcnow().strftime('%Y%m%dT%H%M%SZ')
        out_path = os.path.join(OUT_DIR, f'evidence_tpv30_scec_{ts}.json')
        with open(out_path, 'w') as f:
            json.dump(dict(provenance=prov, submission_res=args.submission_res,
                            submission_dir=submission_dir,
                            reference_file=REFERENCE_FILE,
                            case_params=p, same_physics_verdict=same_physics,
                            **res), f, indent=2, default=str)
        print(f'\nWrote {out_path}')
    return 0


if __name__ == '__main__':
    sys.exit(main())


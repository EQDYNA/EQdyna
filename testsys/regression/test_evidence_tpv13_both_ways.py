#! /usr/bin/env python3
"""
Both-ways demonstration for testsys/parity/evidence_tpv13_scec_comparison.py
(victor-reyes audit, PR #147 (tpv12) and its follow-on tpv13 port, rule 14a): a gate that can only ever print
SUCCESS is not a gate. This drives the real script as a subprocess (never
re-implements its comparison) against:

  1. the real, committed test.reference.results/test.tpv13/ data, compared
     against a SYNTHETIC Barall-format archive built from that same
     committed data (see _build_synthetic_archive) -- CI has no access to
     ~/shared_dataset, so the real archive is never required here; the real
     archive comparison stays the committed evidence script, run by hand
     -> expect exit 0 (every bound met)
  2. the SAME data, with one station file's v-vel column scaled 20x in a
     throwaway tempdir copy (frt.canonical.txt copied unmodified)
     -> expect nonzero exit (peak-ratio bound FAILED)
  3. the same tempdir, with frt.canonical.txt removed entirely
     -> expect nonzero exit, and the message naming the missing file (rule 1:
        raise, not a silent skip)
  4. the synthetic archive's cplot rupture times collapsed to "never
     ruptured" for most non-barrier nodes, stations unmodified
     -> expect nonzero exit, with the FRACTION bound specifically FAILing
        while the station peak-ratio bounds still PASS (a total-collapse
        scenario)
  5. OUR OWN frt.canonical.txt with ~45% of its ruptured non-barrier nodes
     un-ruptured (a PARTIAL STALL, our fraction drops to ~0.44, not a total
     collapse down near 0), synthetic archive unmodified (its own fraction
     stays 0.8051, the real reference value)
     -> expect nonzero exit on the FRACTION_BOUND = 0.28 bound specifically
        (finding 5, PR #155 audit: the original both-ways test only
        exercised a TOTAL collapse of the reference side; a realistic
        partial stall on OUR side, one a human reviewer might plausibly
        mistake for "mostly working" (more than half the fault still
        ruptures), must also fail under the tightened, measured bound.
        Note: this fixture's delta is evaluated against the SYNTHETIC
        archive's fixed 0.8051, not Barall's true 0.9834, so the stall
        fraction needed here to cross FRACTION_BOUND=0.28 is deeper than
        the ~0.60 a run against the real archive would need -- see the
        evidence script's own docstring for the real-archive numbers.)

Case 1 proves the script does not always fail; cases 2-4 prove it does not
always pass, on three different bounds. Together they are the "can come out
both ways" evidence rule 14a requires -- neither a hand-run demo nor prose in
a commit message.

Cheap (rule 9): copies ~100 KB of already-committed station/frt text, no
solver run, no network.
"""
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SCRIPT = os.path.join(REPO_ROOT, 'testsys', 'parity',
                      'evidence_tpv13_scec_comparison.py')
REF_DIR = os.path.join(REPO_ROOT, 'test.reference.results', 'test.tpv13')
DIP_DEG = 60.0  # same TPV13 dip evidence_tpv13_scec_comparison.py uses

STATIONS = [
    ('body010st000dp000', 'body010st000dp000.txt'),
    ('body030st120dp000', 'body030st120dp000.txt'),
]


def _build_synthetic_archive(dest_dir, collapse_rupture=False):
    """A real-shaped Barall-format archive (cplot + the two off-fault body
    station files this script needs), built from this repo's OWN committed
    test.tpv13 reference -- no ~/shared_dataset dependency (finding 1, PR
    #147 audit). Rows match frt.canonical.txt's own (x, z, rupture-time)
    exactly (converted to Barall's (x, down-dip-distance, t) column
    convention), so comparing it against that same reference is, by
    construction, a zero-delta/ratio-1.0 case: a real shape, not a toy grid,
    without needing the real independent archive to prove the script's exit
    codes work both ways.

    collapse_rupture=True sets every non-barrier node's rupture time past
    the gate term, i.e. "never ruptured within the window" -- this is what
    drives the FRACTION bound (not the peak-ratio bound) to FAIL, for the
    both-ways case that targets that bound specifically.
    """
    os.makedirs(dest_dir, exist_ok=True)
    frt = np.loadtxt(os.path.join(REF_DIR, 'frt.canonical.txt'))
    x, z, t = frt[:, 0], frt[:, 2], frt[:, 3]
    downdip = np.abs(z) / math.sin(math.radians(DIP_DEG))
    if collapse_rupture:
        xmin, xmax, dmax = x.min(), x.max(), downdip.max()
        barrier = (np.isclose(x, xmin, atol=1.0) |
                   np.isclose(x, xmax, atol=1.0) |
                   np.isclose(downdip, dmax, atol=1.0))
        t = np.where(barrier, t, 99999.0)
    with open(os.path.join(dest_dir, 'cplot'), 'w') as f:
        f.write('# problem = TPV13 (synthetic, built from this repo\'s own '
                'committed frt.canonical.txt -- see test_evidence_tpv13_'
                'both_ways.py)\n')
        f.write('# Column #1 = horizontal coordinate, distance along '
                'strike (m)\n')
        f.write('# Column #2 = vertical coordinate, distance down-dip (m)\n')
        f.write('# Column #3 = rupture time (s)\n')
        for xi, di, ti in zip(x, downdip, t):
            f.write('%.6e %.6e %.6e\n' % (xi, di, ti))
    for archive_name, our_name in STATIONS:
        shutil.copy(os.path.join(REF_DIR, 'stations', our_name),
                    os.path.join(dest_dir, archive_name))


def _build_partial_stall_frt(dest_path, stall_fraction=0.45):
    """A copy of the real committed frt.canonical.txt with STALL_FRACTION of
    its ruptured, non-barrier nodes un-ruptured (rupture time set past the
    gate term) -- a PARTIAL stall (finding 5, PR #155 audit), distinct from
    the existing total-collapse scenario in _build_synthetic_archive. Uses
    the same barrier convention as the real script (exclude the two
    along-strike edges and the deepest down-dip row)."""
    frt = np.loadtxt(os.path.join(REF_DIR, 'frt.canonical.txt'))
    x, z, t = frt[:, 0], frt[:, 2], frt[:, 3]
    downdip = np.abs(z) / math.sin(math.radians(DIP_DEG))
    xmin, xmax, dmax = x.min(), x.max(), downdip.max()
    barrier = (np.isclose(x, xmin, atol=1.0) | np.isclose(x, xmax, atol=1.0) |
               np.isclose(downdip, dmax, atol=1.0))
    ruptured = (~barrier) & (t < 999.0)
    idx = np.flatnonzero(ruptured)
    # Deterministic selection (every Nth ruptured node), not random -- a
    # both-ways test must be reproducible.
    n_stall = int(round(stall_fraction * idx.size))
    stall_idx = idx[np.linspace(0, idx.size - 1, n_stall, dtype=int)]
    out = frt.copy()
    out[stall_idx, 3] = 99999.0
    np.savetxt(dest_path, out, fmt='%.6e')


def _run(run_dir, gate_resolution=True, archive_dir=None):
    cmd = [sys.executable, SCRIPT, '--run-dir', run_dir, '--no-json']
    if gate_resolution:
        cmd.insert(-1, '--gate-resolution')
    env = dict(os.environ)
    if archive_dir is not None:
        env['TPV13_BARALL_ARCHIVE_DIR'] = archive_dir
    p = subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True,
                        env=env)
    return p.returncode, p.stdout + p.stderr


def _scale_v_vel(path, factor):
    lines = open(path).readlines()
    out = []
    for line in lines:
        cols = line.strip().split()
        is_data = bool(cols)
        if is_data:
            try:
                float(cols[0])
            except ValueError:
                is_data = False
        if not is_data:
            out.append(line)
            continue
        cols[4] = repr(float(cols[4]) * factor)
        out.append('  '.join(cols) + '\n')
    open(path, 'w').writelines(out)


def main():
    if not os.path.isdir(REF_DIR):
        raise SystemExit('%s: missing -- this test requires the committed '
                         'test.tpv13 reference tree' % REF_DIR)

    fails = []
    tmp = tempfile.mkdtemp(prefix='tpv13_both_ways_')
    try:
        archive = os.path.join(tmp, 'synthetic_barall_archive')
        _build_synthetic_archive(archive)

        # 1. real reference data vs a synthetic-but-real-shaped archive built
        #    from that same reference -> PASS (exit 0). No ~/shared_dataset
        #    dependency (finding 1, PR #147 (tpv12) and its follow-on tpv13 port audit): CI never fetches that
        #    data, so the real archive must not be required for this test to
        #    run at all.
        rc, out = _run(os.path.join(REF_DIR, 'stations'), archive_dir=archive)
        ok = (rc == 0)
        print('  real reference data                : rc=%d  %s' % (rc, 'ok' if ok else 'FAIL (expected 0)'))
        if not ok:
            fails.append('real reference data did not exit 0:\n' + out)

        # 2. perturbed v-vel -> FAIL (nonzero exit, peak-ratio bound)
        stations_tmp = os.path.join(tmp, 'stations')
        shutil.copytree(os.path.join(REF_DIR, 'stations'), stations_tmp)
        shutil.copy(os.path.join(REF_DIR, 'frt.canonical.txt'), tmp)
        _scale_v_vel(os.path.join(stations_tmp, 'body010st000dp000.txt'), 20.0)
        rc, out = _run(tmp, archive_dir=archive)
        ok = (rc != 0) and ('FAIL' in out)
        print('  v-vel scaled 20x                   : rc=%d  %s' % (rc, 'ok' if ok else 'FAIL (expected nonzero + FAIL verdict)'))
        if not ok:
            fails.append('perturbed v-vel did not fail as expected:\n' + out)

        # 3. missing frt.canonical.txt -> raises (nonzero exit, named file)
        os.remove(os.path.join(tmp, 'frt.canonical.txt'))
        rc, out = _run(tmp, archive_dir=archive)
        ok = (rc != 0) and ('frt.canonical.txt' in out) and ('missing' in out)
        print('  frt.canonical.txt removed          : rc=%d  %s' % (rc, 'ok' if ok else 'FAIL (expected a named-file raise)'))
        if not ok:
            fails.append('missing frt.canonical.txt did not raise as expected:\n' + out)

        # 4. fraction bound specifically FAILs: synthetic archive's rupture
        #    collapsed to "unruptured" everywhere but the (unbreakable)
        #    barrier nodes, real stations unmodified -> nonzero exit, with
        #    the FRACTION_BOUND line FAILing and the station peak-ratio
        #    bounds still PASSing (finding 5, PR #147 (tpv12) and its follow-on tpv13 port audit: a case that
        #    exercises THIS bound, not just the peak-ratio one above).
        archive_frac_fail = os.path.join(tmp, 'synthetic_barall_archive_fracfail')
        _build_synthetic_archive(archive_frac_fail, collapse_rupture=True)
        rc, out = _run(os.path.join(REF_DIR, 'stations'),
                        archive_dir=archive_frac_fail)
        frac_line = next((l for l in out.splitlines()
                           if 'bound FRACTION_BOUND' in l), '')
        ok = (rc != 0) and frac_line.rstrip().endswith('FAIL') and \
            'all station peak-ratio bounds: PASS' in out
        print('  fraction bound collapsed           : rc=%d  %s' % (rc, 'ok' if ok else 'FAIL (expected FRACTION_BOUND FAIL, station bounds PASS)'))
        if not ok:
            fails.append('collapsed-rupture fixture did not fail the '
                          'fraction bound (and only that bound) as '
                          'expected:\n' + out)

        # 5. PARTIAL stall on OUR side (our fraction drops to ~0.44, not a
        #    total collapse down near 0): archive unmodified (its own real
        #    fraction, 0.8051), 45% of our ruptured non-barrier nodes
        #    un-ruptured -> fraction ours drops from 0.8051 to ~0.44,
        #    |delta| vs the archive's 0.8051 lands well past
        #    FRACTION_BOUND = 0.28 (finding 5, PR #155 audit). Checking
        #    frac_ours > 0.3 (rather than near 0) is what tells this apart
        #    from scenario 4's TOTAL collapse.
        stall_dir = os.path.join(tmp, 'stall')
        os.makedirs(stall_dir, exist_ok=True)
        shutil.copytree(os.path.join(REF_DIR, 'stations'),
                         os.path.join(stall_dir, 'stations'))
        _build_partial_stall_frt(os.path.join(stall_dir, 'frt.canonical.txt'))
        rc, out = _run(stall_dir, archive_dir=archive)
        frac_line = next((l for l in out.splitlines()
                           if 'bound FRACTION_BOUND' in l), '')
        frac_ours_line = next((l for l in out.splitlines()
                                if l.strip().startswith('ours  :')), '')
        m = re.search(r'\(([0-9.]+)\)', frac_ours_line)
        frac_ours = float(m.group(1)) if m else None
        ok = (rc != 0) and frac_line.rstrip().endswith('FAIL') and \
            frac_ours is not None and 0.30 < frac_ours < 0.55
        print('  partial stall (fraction ~0.44)     : rc=%d  %s' % (rc, 'ok' if ok else 'FAIL (expected FRACTION_BOUND FAIL at a partial, not total, fraction)'))
        if not ok:
            fails.append('partial-stall fixture did not fail the fraction '
                          'bound at a partial (not total) fraction as '
                          'expected:\n' + out)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    if fails:
        print('\nFAIL test_evidence_tpv13_both_ways (%d)' % len(fails))
        for f in fails:
            print(' -', f)
        return 1
    print('\nSUCCESS test_evidence_tpv13_both_ways: PASS on real data vs a '
          'synthetic real-shaped archive, FAIL on a perturbed station file, '
          'raise on a missing frt.canonical.txt, FAIL on the fraction bound '
          'alone with a total-collapse fixture, FAIL on the fraction bound '
          'alone with a partial-stall (fraction ~0.44, not total) fixture')
    return 0


if __name__ == '__main__':
    sys.exit(main())

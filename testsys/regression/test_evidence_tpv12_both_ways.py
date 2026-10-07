#! /usr/bin/env python3
"""
Both-ways demonstration for testsys/parity/evidence_tpv12_scec_comparison.py
(victor-reyes audit, PR #147, rule 14a): a gate that can only ever print
SUCCESS is not a gate. This drives the real script as a subprocess (never
re-implements its comparison) against:

  1. the real, committed test.reference.results/test.tpv12/ data
     -> expect exit 0 (every bound met)
  2. the SAME data, with one station file's v-vel column scaled 20x in a
     throwaway tempdir copy (frt.canonical.txt copied unmodified)
     -> expect nonzero exit (peak-ratio bound FAILED)
  3. the same tempdir, with frt.canonical.txt removed entirely
     -> expect nonzero exit, and the message naming the missing file (rule 1:
        raise, not a silent skip)

Case 1 proves the script does not always fail; cases 2-3 prove it does not
always pass. Together they are the "can come out both ways" evidence rule
14a requires -- neither a hand-run demo nor prose in a commit message.

Cheap (rule 9): copies ~100 KB of already-committed station text, no run,
no network.
"""
import os
import shutil
import subprocess
import sys
import tempfile

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SCRIPT = os.path.join(REPO_ROOT, 'testsys', 'parity',
                      'evidence_tpv12_scec_comparison.py')
REF_DIR = os.path.join(REPO_ROOT, 'test.reference.results', 'test.tpv12')


def _run(run_dir, gate_resolution=True):
    cmd = [sys.executable, SCRIPT, '--run-dir', run_dir, '--no-json']
    if gate_resolution:
        cmd.insert(-1, '--gate-resolution')
    p = subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True)
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
                         'test.tpv12 reference tree' % REF_DIR)

    fails = []

    # 1. real data -> PASS (exit 0)
    rc, out = _run(os.path.join(REF_DIR, 'stations'))
    ok = (rc == 0)
    print('  real reference data                : rc=%d  %s' % (rc, 'ok' if ok else 'FAIL (expected 0)'))
    if not ok:
        fails.append('real reference data did not exit 0:\n' + out)

    tmp = tempfile.mkdtemp(prefix='tpv12_both_ways_')
    try:
        # 2. perturbed v-vel -> FAIL (nonzero exit)
        stations_tmp = os.path.join(tmp, 'stations')
        shutil.copytree(os.path.join(REF_DIR, 'stations'), stations_tmp)
        shutil.copy(os.path.join(REF_DIR, 'frt.canonical.txt'), tmp)
        _scale_v_vel(os.path.join(stations_tmp, 'body010st000dp000.txt'), 20.0)
        rc, out = _run(tmp)
        ok = (rc != 0) and ('FAIL' in out)
        print('  v-vel scaled 20x                   : rc=%d  %s' % (rc, 'ok' if ok else 'FAIL (expected nonzero + FAIL verdict)'))
        if not ok:
            fails.append('perturbed v-vel did not fail as expected:\n' + out)

        # 3. missing frt.canonical.txt -> raises (nonzero exit, named file)
        os.remove(os.path.join(tmp, 'frt.canonical.txt'))
        rc, out = _run(tmp)
        ok = (rc != 0) and ('frt.canonical.txt' in out) and ('missing' in out)
        print('  frt.canonical.txt removed          : rc=%d  %s' % (rc, 'ok' if ok else 'FAIL (expected a named-file raise)'))
        if not ok:
            fails.append('missing frt.canonical.txt did not raise as expected:\n' + out)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    if fails:
        print('\nFAIL test_evidence_tpv12_both_ways (%d)' % len(fails))
        for f in fails:
            print(' -', f)
        return 1
    print('\nSUCCESS test_evidence_tpv12_both_ways: PASS on real data, FAIL on '
          'a perturbed station file, raise on a missing frt.canonical.txt')
    return 0


if __name__ == '__main__':
    sys.exit(main())

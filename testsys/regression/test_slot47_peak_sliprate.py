#! /usr/bin/env python3
"""
Regression guard (PROJECT_RULES.md rule 10) for the slot-47 slip-rate bug,
pathway_forward.md item 17(b) (owner decision 2026-10-02, PR #76).

THE BUG: FRIC_SLOT_PEAK_SLIPRATE (fric slot 47, globalvar.f90:57) was written
ONLY inside solveRSF (src/fortran/faulting.f90:275; python faulting.py, the
`fric = B.setat(xp, fric, (slice(None), gv.PEAK_SLIPRATE), v_trial)` line in
solveRSF), which only ever runs for friclaw>=3 (rate-and-state). For friclaw
1/2 (slip-weakening/time-weakening -- solveSWTW, never touches this slot),
the slot stayed at its restart-init value (readInputFiles.py:480,
netcdf_io.f90:110/185 -- typically 0 or the creep floor) for the entire run.
That zero then propagated to:
  - frt.txt's slip-rate output column (library_output.f90:349, "47: final
    slip rate" -- column index 10, 0-based, in the 22-column canonical form)
  - src_evol's final slip rate (library_output.f90:506)
  - the restart netCDF slip rate

MEASURED, before the fix (test.reference.results/test.tpv8/frt.canonical.txt,
committed under the old code): 830 of 1891 fault nodes ruptured (fnft in
[0.087, 4.986] s, well inside this case's 5 s term) and EVERY ONE of them
reports peak-sliprate column == 0.0 -- physically impossible for a node that
ruptured.

THE FIX (this commit): getNsdSlipSliprateTraction
(src/fortran/faulting.f90:106-115; src/python/eqdyna/faulting.py, right after
srMag/CUM_SLIP are written) now writes FRIC_SLOT_PEAK_SLIPRATE from
nsdSliprateVector(4)/srMag (the 3-component slip-rate magnitude) for EVERY
friction law, every step. solveRSF's own write (friclaw>=3 only) still runs
immediately after and overwrites it with v_trial there, so friclaw>=3 output
is provably unchanged by this fix (verified separately, as its own commit,
against drv.a6/tpv104/tpv1053d's committed references -- byte-identical).
This test only needs to prove the friclaw<3 (slip-weakening) case, which is
the one this bug actually silently broke.

Fixture: test.tpv8 (friclaw=1, the project's own basic slip-weakening gate
case), forced serial (npx=npy=npz=1) with a short par.term=1.0 s -- test.tpv8's
own rupture times start at 0.087 s (measured above), so 1 s already covers
many ruptured nodes while running in a few seconds, not this case's full 5 s
gate term.

Exits non-zero on any failure.
"""
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpirun_capture  # noqa: E402

if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from testsys import frt_canonical  # noqa: E402

MPIRUN = os.environ.get('EQDYNA_MPIRUN', 'mpirun')
CASE_NAME = 'test.tpv8'
SHORT_TERM_S = 1.0  # covers rupture times down to 0.087s (measured); full gate term is 5s.
PEAK_SLIPRATE_COL = 10  # 0-based; "47: final slip rate" in library_output.f90:349's write order.
FNFT_COL = 3            # 0-based; rupture time.
FNFT_UNRUPTURED_ABOVE = 5000.0


def _env():
    env = dict(os.environ)
    env['EQDYNAROOT'] = ROOT
    env['PATH'] = os.pathsep.join([os.path.join(ROOT, 'bin'),
                                   os.path.join(ROOT, 'scripts'),
                                   env.get('PATH', '')])
    return env


def _require_binary():
    exe = os.path.join(ROOT, 'bin', 'eqdyna')
    if not os.path.isfile(exe):
        raise AssertionError(
            'bin/eqdyna is missing -- build it first (./install-eqdyna.sh -m <platform>) '
            'before running this test.')
    return exe


def _build_and_run_serial(tmp):
    case_dir = os.path.join(tmp, 'case_serial')
    env = _env()
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'create.newcase'),
                        case_dir, CASE_NAME], env=env, capture_output=True, text=True)
    if r.returncode != 0:
        raise AssertionError('create.newcase failed (%d):\n%s' % (r.returncode, r.stdout + r.stderr))
    params = os.path.join(case_dir, 'user_defined_params.py')
    with open(params) as f:
        text = f.read().rstrip('\n')
    with open(params, 'w') as f:
        f.write(text + '\n\n# forced serial + shortened term by this test (rule 9: cheap)\n'
                       'par.nx = 1\npar.ny = 1\npar.nz = 1\n'
                       'par.term = %r\n' % SHORT_TERM_S)
    r = subprocess.run([sys.executable, 'case.setup'], cwd=case_dir, env=env,
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise AssertionError('case.setup failed (%d):\n%s' % (r.returncode, r.stdout + r.stderr))

    exe = _require_binary()
    rc, out = mpirun_capture.run_rank_files(MPIRUN, exe, case_dir, np=1, env=env,
                                            timeout=300)
    if rc != 0:
        raise AssertionError('eqdyna exited %d for serial test.tpv8:\n%s' % (rc, out[-3000:]))
    return case_dir


def check_slipweakening_frt_sliprate_column_nonzero_after_rupture():
    """friclaw=1 (test.tpv8): every node that ruptured (fnft < sentinel) must
    report a nonzero peak-slip-rate column. RED on the pre-fix code: measured
    830/1891 ruptured nodes, all reporting exactly 0.0 in this column (see
    module docstring)."""
    with tempfile.TemporaryDirectory() as tmp:
        case_dir = _build_and_run_serial(tmp)
        rows = frt_canonical.canonical_from_case(case_dir)
    fnft = rows[:, FNFT_COL]
    ruptured_mask = fnft < FNFT_UNRUPTURED_ABOVE
    n_ruptured = int(ruptured_mask.sum())
    if n_ruptured == 0:
        raise AssertionError(
            'serial test.tpv8 at term=%.1fs: NO node ruptured (all fnft >= sentinel) -- '
            'fixture does not exercise the bug at all, strengthen it (nucleation timing?)'
            % SHORT_TERM_S)
    sliprate = rows[ruptured_mask, PEAK_SLIPRATE_COL]
    n_zero = int((sliprate == 0.0).sum())
    if n_zero > 0:
        raise AssertionError(
            'serial test.tpv8 (friclaw=1, slip-weakening): %d of %d ruptured fault nodes '
            'report EXACTLY 0.0 in the peak-slip-rate frt column (col %d) -- slot-47 bug '
            '(FRIC_SLOT_PEAK_SLIPRATE written only in solveRSF, friclaw>=3) is back. '
            'min/max sliprate among ruptured: %r/%r'
            % (n_zero, n_ruptured, PEAK_SLIPRATE_COL, float(sliprate.min()), float(sliprate.max())))
    print('  PASS  serial test.tpv8 (friclaw=1): %d/%d ruptured fault nodes all report '
          'nonzero peak slip rate (min %.6e, max %.6e m/s)'
          % (n_ruptured, n_ruptured, float(sliprate.min()), float(sliprate.max())))


def main():
    print('Regression guard: slot-47 (FRIC_SLOT_PEAK_SLIPRATE) written for every friclaw (item 17b)')
    failures = []
    for c in (check_slipweakening_frt_sliprate_column_nonzero_after_rupture,):
        try:
            c()
        except AssertionError as e:
            failures.append('%s: %s' % (c.__name__, e))
            print('  FAIL  %s: %s' % (c.__name__, e))
    if failures:
        print('\nFAIL test_slot47_peak_sliprate (%d check(s))' % len(failures))
        return 1
    print('\nSUCCESS test_slot47_peak_sliprate')
    return 0


if __name__ == '__main__':
    sys.exit(main())

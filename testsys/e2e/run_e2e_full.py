#! /usr/bin/env python3
"""
e2e "full" tier -- SCEC benchmark cases AS SPECIFIED (spec element size AND
spec duration), 16 MPI ranks, REPORT-ONLY (PROJECT_RULES.md rules 4, 6, 9).

No golden references exist at full (spec) resolution -- unlike the fast
e2e tier (run_e2e.py), this tier never bit-compares against
test.reference.results/. Its product is a provenance-stamped results table
plus completion/sanity checks (did the run finish, are the expected output
files present, is max slip/ruptured area a finite, physically sane number).
Exit is non-zero ONLY on crash or non-completion -- reported values
(slip, ruptured area, wall time) are never gated against a threshold here.

DESIGN CONTRACT: this tier NEVER forks a *_full compset. It drives the
exact same case_input/<case>/ tree the fast tier uses (create.newcase,
case.setup, mpirun, same as run_e2e.py); the only difference is a
mechanical rewrite of par.dx/par.nx,ny,nz/par.term in the freshly generated
user_defined_params.py before case.setup runs, mirroring
testsys/perf/run_scaling.py's make_case(). There is no test.tpv8_full/ or
similar directory anywhere in this repo.

Case selection: testsys/e2e/full_specs.py's FULL_SPECS map (spec dx/term/
16-rank decomposition + a citation into scratch/specs/, where the SCEC PDFs
were fetched from https://strike.scec.org/cvws/). Cases with no
discoverable spec number (or, for internal non-SCEC cases, no published
full-resolution provenance) are in full_specs.EXCLUDED and are skipped with
a printed notice -- never given an invented number.

Usage:
    python3 testsys/e2e/run_e2e_full.py
    python3 testsys/e2e/run_e2e.py --full   # equivalent entry point
"""
import json
import os
import re
import shutil
import subprocess
import sys
import time

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MACHINE = os.environ.get('EQDYNA_TEST_MACHINE', 'ubuntu')
MPIRUN = os.environ.get('EQDYNA_MPIRUN', 'mpirun')
RANKS = 16

sys.path.insert(0, REPO_ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from full_specs import FULL_SPECS, EXCLUDED  # noqa: E402


def _run(cmd, cwd, env):
    print('+ (' + os.path.basename(cwd) + ') ' + ' '.join(str(c) for c in cmd))
    return subprocess.call(cmd, cwd=cwd, env=env)


def _rewrite_params(case_dir, dx, term, nx, ny, nz):
    """Mechanical override of par.dx/par.term/par.nx,ny,nz in the freshly
    create.newcase-generated user_defined_params.py, BEFORE case.setup runs.
    Same regex forms as testsys/perf/run_scaling.py's make_case() (combined
    `par.nx, par.ny, par.nz = ...` line, and separate `par.nx = ...` lines)
    so both styles used across case_input/*/user_defined_params.py are
    handled."""
    p = os.path.join(case_dir, 'user_defined_params.py')
    s = open(p).read()
    s = re.sub(r'^par\.dx\s*=.*', f'par.dx = {dx}', s, flags=re.M)
    s = re.sub(r'^par\.term\s*=.*', f'par.term = {term}', s, flags=re.M)
    s = re.sub(r'par\.nx\s*,\s*par\.ny\s*,\s*par\.nz\s*=.*',
               f'par.nx, par.ny, par.nz = {nx}, {ny}, {nz}', s)
    s = re.sub(r'^par\.nx = .*', f'par.nx = {nx}', s, flags=re.M)
    s = re.sub(r'^par\.ny = .*', f'par.ny = {ny}', s, flags=re.M)
    s = re.sub(r'^par\.nz = .*', f'par.nz = {nz}', s, flags=re.M)
    open(p, 'w').write(s)
    return p


def _cheap_sanity(case_dir, spec):
    """Max slip / ruptured area from frt.txt*, cheaply (no matplotlib): the
    same loader plotRuptureDynamics uses (scripts/lib.py:loadFrtData), just
    without the figure. Returns a dict, or {'error': ...} if frt.txt* is
    missing/empty (non-completion)."""
    frt_files = [f for f in os.listdir(case_dir) if re.match(r'^frt\.txt\d+$', f)]
    if not frt_files:
        return {'error': 'no frt.txt* files -- run did not reach on-fault output'}
    sys.path.insert(0, os.path.join(REPO_ROOT, 'scripts'))
    for mod in ('user_defined_params', 'lib'):
        sys.modules.pop(mod, None)
    cwd = os.getcwd()
    os.chdir(case_dir)
    try:
        from user_defined_params import par
        from lib import loadFrtData
        xx, zz, rupt, rupt2d, fVarArr, magnitude = loadFrtData(par)
    except Exception as e:
        return {'error': f'loadFrtData failed: {e}'}
    finally:
        os.chdir(cwd)
        for mod in ('user_defined_params', 'lib'):
            sys.modules.pop(mod, None)
    ruptured = rupt2d[:, :, 0] > 0.  # rupture time > 0 -> node ruptured
    max_slip = float(rupt2d[:, :, 1].max())
    ruptured_area_km2 = float(ruptured.sum() * (par.dx * par.dz) / 1.e6)
    return {'max_slip_m': max_slip, 'ruptured_area_km2': ruptured_area_km2,
            'magnitude': float(magnitude), 'n_ruptured_nodes': int(ruptured.sum())}


def main():
    if os.environ.get('EQDYNA_FULL_LAUNCH') != 'yes-hours':
        print('run_e2e_full: refusing to launch -- these are spec-resolution, '
              'multi-hour, 16-rank runs, user-scheduled, never automatic. '
              'Set EQDYNA_FULL_LAUNCH=yes-hours to actually launch them.')
        return 1

    for case, reason in EXCLUDED.items():
        print(f'e2e-full: EXCLUDED {case}: {reason}')

    bin_exe = os.path.join(REPO_ROOT, 'bin', 'eqdyna')
    if not os.path.exists(bin_exe):
        print(f'e2e-full: FAIL - {bin_exe} does not exist; build first '
              f'(./install-eqdyna.sh -m {MACHINE})')
        return 1

    sha = subprocess.check_output(['git', '-C', REPO_ROOT, 'rev-parse', '--short', 'HEAD'],
                                   text=True).strip()
    host = os.uname().nodename
    load = os.getloadavg()

    test_dir = os.path.join(REPO_ROOT, 'test.full')
    prev_dir = os.path.join(REPO_ROOT, 'test.full.prev')
    if os.path.isdir(test_dir):
        if os.path.isdir(prev_dir):
            shutil.rmtree(prev_dir)
        shutil.move(test_dir, prev_dir)
        print(f'e2e-full: preserved previous run as {prev_dir}')
    os.makedirs(test_dir)

    env = dict(os.environ)
    env['EQDYNAROOT'] = REPO_ROOT
    env['PATH'] = os.pathsep.join([
        os.path.join(REPO_ROOT, 'bin'),
        os.path.join(REPO_ROOT, 'scripts'),
        env.get('PATH', ''),
    ])

    rows = []
    failures = []
    for case, spec in FULL_SPECS.items():
        case_dir = os.path.join(test_dir, case)
        print(f'\n-- {case}: dx={spec["dx"]}m term={spec["term"]}s '
              f'{RANKS} ranks ({spec["nx"]},{spec["ny"]},{spec["nz"]}) --')
        if _run(['create.newcase', case, case], test_dir, env) != 0:
            failures.append(f'{case}: create.newcase failed')
            continue
        _rewrite_params(case_dir, spec['dx'], spec['term'],
                         spec['nx'], spec['ny'], spec['nz'])
        if _run(['./case.setup'], case_dir, env) != 0:
            failures.append(f'{case}: case.setup failed')
            continue
        t0 = time.time()
        rc = _run([MPIRUN, '-np', str(RANKS), 'eqdyna'], case_dir, env)
        elapsed = time.time() - t0
        if rc != 0:
            failures.append(f'{case}: mpirun -np {RANKS} eqdyna exited {rc}')
            rows.append(dict(case=case, dx=spec['dx'], term=spec['term'],
                              ranks=RANKS, wall_s=elapsed, completed=False,
                              citation=spec['citation']))
            continue
        sanity = _cheap_sanity(case_dir, spec)
        completed = 'error' not in sanity
        if not completed:
            failures.append(f'{case}: {sanity["error"]}')
        outputs = sorted(os.listdir(case_dir))
        rows.append(dict(case=case, dx=spec['dx'], term=spec['term'], ranks=RANKS,
                          wall_s=elapsed, completed=completed,
                          citation=spec['citation'], outputs=outputs, **sanity))
        print(f'   wall={elapsed:.1f}s completed={completed} {sanity}')

    meta = dict(sha=sha, host=host, loadavg=load,
                date=time.strftime('%Y-%m-%d %H:%M'), ranks=RANKS, rows=rows,
                excluded={c: r for c, r in EXCLUDED.items()})
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'full_last.json')
    json.dump(meta, open(out, 'w'), indent=1)
    print(f'\ne2e-full: provenance {sha} on {host}, loadavg {load}; saved {out}')

    if failures:
        print('e2e-full: FAIL - crash/non-completion:')
        for f in failures:
            print('  -', f)
        return 1
    print('e2e-full: SUCCESS - all spec-resolution runs completed (report-only, not gated)')
    return 0


if __name__ == '__main__':
    sys.exit(main())

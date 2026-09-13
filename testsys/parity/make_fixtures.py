#! /usr/bin/env python3
"""
Builds `eqdyna-pydump` (src/makefile's PYDUMP=1 target) and generates the
tpv8 serial parity fixtures the `parity` testsys tier compares the Python
ports against: a fresh Fortran `frt.txt` (the golden oracle) plus the
`pydump_*` static-state and step1-5 checkpoint dumps `python/eqdyna/port.py`
loads directly.

Regenerable, not a static blob to commit: this script IS the fixture
provenance (rule 4, "only fresh runs are evidence"). Run it whenever
src/*.f90 changes; the `parity` tier gates against whatever fixtures are
on disk, so a stale fixture is a self-inflicted false pass -- re-run this
first if in doubt.

Usage: python3 testsys/parity/make_fixtures.py [--case test.tpv8] [--term 5.0]
"""
import argparse
import os
import shutil
import subprocess
import sys

TESTSYS = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(TESTSYS))
FIXTURE_ROOT = os.path.join(TESTSYS, 'fixtures')


def sh(cmd, cwd=None, env=None):
    print('+', ' '.join(cmd), (f'(cwd={cwd})' if cwd else ''))
    r = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-4000:]); print(r.stderr[-4000:])
        raise SystemExit(f'FAIL: {" ".join(cmd)} exited {r.returncode}')
    return r


def build_pydump_binary():
    src = os.path.join(REPO_ROOT, 'src')
    env = dict(os.environ, MACHINE=os.environ.get('MACHINE', 'ubuntu'))
    # Clean rebuild of BOTH targets so a stale default `eqdyna.o` set left
    # over from a previous build can't silently get reused by either.
    subprocess.run(['bash', '-c', 'rm -f *.o eqdyna eqdyna-pydump'], cwd=src)
    sh(['make', f'MACHINE={env["MACHINE"]}'], cwd=src, env=env)
    default_bin = os.path.join(src, 'eqdyna')
    default_copy = os.path.join(FIXTURE_ROOT, 'eqdyna-default-for-neutrality-check')
    os.makedirs(FIXTURE_ROOT, exist_ok=True)
    shutil.copy(default_bin, default_copy)
    sh(['make', f'MACHINE={env["MACHINE"]}', 'eqdyna-pydump'], cwd=src, env=env)
    pydump_bin = os.path.join(src, 'eqdyna-pydump')
    if not os.path.exists(pydump_bin):
        raise SystemExit('FAIL: src/eqdyna-pydump was not produced')
    return pydump_bin, default_copy


def make_case(case_name, case_dir, term):
    env = dict(os.environ, EQDYNAROOT=REPO_ROOT,
               PYTHONPATH=os.path.join(REPO_ROOT, 'scripts'))
    if os.path.isdir(case_dir):
        shutil.rmtree(case_dir)
    sh([sys.executable, os.path.join(REPO_ROOT, 'scripts', 'create.newcase'),
        case_dir, case_name], env=env)
    params_path = os.path.join(case_dir, 'user_defined_params.py')
    with open(params_path, 'a') as f:
        f.write('\npar.nx = 1\npar.ny = 1\npar.nz = 1\n')
        if term is not None:
            f.write(f'par.term = {term}\n')
    r = subprocess.run([sys.executable, 'case.setup'], cwd=case_dir, env=env,
                        capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-4000:]); print(r.stderr[-4000:])
        raise SystemExit('FAIL: case.setup')
    return case_dir


def run_case(binary, case_dir):
    dst = os.path.join(case_dir, os.path.basename(binary))
    shutil.copy(binary, dst)
    r = subprocess.run(['mpirun', '-np', '1', '-wdir', case_dir, dst],
                        capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-4000:]); print(r.stderr[-4000:])
        raise SystemExit(f'FAIL: {dst} exited {r.returncode} on {case_dir}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', default='test.tpv8')
    ap.add_argument('--term', type=float, default=None,
                     help='override par.term (smaller -> fewer steps, cheaper fixture)')
    args = ap.parse_args()

    pydump_bin, default_bin = build_pydump_binary()

    fixture_case = os.path.join(FIXTURE_ROOT, args.case.replace('.', '_') + '_serial')
    make_case(args.case, fixture_case, args.term)
    run_case(pydump_bin, fixture_case)

    required = ['frt.txt0', 'pydump_header.txt', 'pydump_meshCoor.txt', 'pydump_conn.txt',
                'pydump_elemgeo.txt', 'pydump_nodeinfo.txt', 'pydump_nodalmass.txt',
                'pydump_fnms.txt', 'pydump_fault.txt', 'pydump_v1.txt', 'pydump_stations.txt',
                'pydump_step5_accel.txt', 'pydump_step5_veldisp.txt', 'pydump_step5_fault.txt']
    missing = [f for f in required if not os.path.exists(os.path.join(fixture_case, f))]
    if missing:
        raise SystemExit(f'FAIL: fixture generation ran but is missing {missing}')

    print(f'\nFixtures written to {fixture_case}')
    print(f'Default (no-op stub) binary saved to {default_bin} for the neutrality check.')
    print('SUCCESS make_fixtures')


if __name__ == '__main__':
    main()

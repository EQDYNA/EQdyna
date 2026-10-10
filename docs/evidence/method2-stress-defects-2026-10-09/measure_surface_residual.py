#! /usr/bin/env python3
"""
Defect C step-1 surface residual acceleration: before (+7.3215 m offset) vs
after (offset removed), serial fortran, gate dx, for the 4 affected cases.

Builds a scratch serial case per name (run_e2e.make_serial_case, rule 1),
runs bin/eqdyna with EQDYNA_DUMP_EQUIL=1, kills it as soon as
equilibriumDump.0000.txt stops growing (dumpNodalAccel writes once at
nt==1; waiting out the rest of the gate term buys nothing for this
measurement), and reports BOTH mean|a_z| and the signed mean a_z over
free-surface (z>-1e-6 m), non-fault (faultTag==0), non-PML (dof<=3) nodes.

"before" is measured by reverting the single call-site line in
src/fortran/meshgen.f90, rebuilding bin/eqdyna, measuring, then restoring
the committed source and rebuilding again -- verified via git diff --stat
showing no change before exiting. Never leaves the tree dirty.
"""
import os
import subprocess
import sys
import tempfile
import time

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))
SCRATCH = os.environ.get('MEASURE_RESIDUAL_SCRATCH') or tempfile.mkdtemp(
    prefix='eqdyna_residual_cases_')
MESHGEN = os.path.join(REPO_ROOT, 'src', 'fortran', 'meshgen.f90')
BIN = os.path.join(REPO_ROOT, 'bin', 'eqdyna')

CASES = ['test.tpv13', 'test.tpv27', 'test.tpv30', 'test.drv.a6']

AFTER_LINE = '                    if (C_elastic == 0) call setPlasticStress(-0.5d0*(zline(iz)+zline(iz-1)), elemCount)\n'
BEFORE_LINE = '                    if (C_elastic == 0) call setPlasticStress(-0.5d0*(zline(iz)+zline(iz-1)) + 7.3215d0, elemCount)\n'


def rebuild():
    env = dict(os.environ, MACHINE='ubuntu')
    r = subprocess.run(['make'], cwd=os.path.join(REPO_ROOT, 'src', 'fortran'), env=env,
                        capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError('make failed:\n%s\n%s' % (r.stdout[-4000:], r.stderr[-4000:]))
    built = os.path.join(REPO_ROOT, 'src', 'fortran', 'eqdyna')
    os.replace(built, BIN)


def set_offset_present(present):
    with open(MESHGEN) as f:
        text = f.read()
    target_from = BEFORE_LINE if present else AFTER_LINE
    target_to = AFTER_LINE if present else BEFORE_LINE
    if target_to not in text:
        if target_from in text:
            return  # already in the requested state
        raise RuntimeError('neither call-site line found in meshgen.f90 -- source has drifted')
    text = text.replace(target_to, target_from)
    with open(MESHGEN, 'w') as f:
        f.write(text)


def build_case(case_name):
    sys.path.insert(0, os.path.join(REPO_ROOT, 'testsys', 'e2e'))
    import run_e2e
    case_dir = os.path.join(SCRATCH, case_name)
    if not os.path.isdir(case_dir):
        os.makedirs(SCRATCH, exist_ok=True)
        run_e2e.make_serial_case(case_name, case_dir, run_e2e.base_env())
    return case_dir


def run_dump(case_dir):
    dump_path = os.path.join(case_dir, 'equilibriumDump.0000.txt')
    if os.path.exists(dump_path):
        os.remove(dump_path)
    env = dict(os.environ, EQDYNA_DUMP_EQUIL='1')
    proc = subprocess.Popen([BIN], cwd=case_dir, env=env,
                             stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    last_size = -1
    stable = 0
    deadline = time.time() + 300
    try:
        while time.time() < deadline:
            time.sleep(1.0)
            if os.path.exists(dump_path):
                size = os.path.getsize(dump_path)
                if size == last_size and size > 0:
                    stable += 1
                    if stable >= 2:
                        break
                else:
                    stable = 0
                last_size = size
    finally:
        proc.kill()
        proc.wait()
    if not os.path.exists(dump_path):
        raise RuntimeError('%s: equilibriumDump.0000.txt never appeared within 300s' % case_dir)
    return dump_path


def stats(dump_path):
    xs, zs, azs, tags, dofs = [], [], [], [], []
    with open(dump_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) != 8:
                continue
            x, y, z, ax, ay, az, tag, dof = parts
            zs.append(float(z))
            azs.append(float(az))
            tags.append(int(tag))
            dofs.append(int(dof))
    import numpy as np
    z = np.array(zs)
    az = np.array(azs)
    tag = np.array(tags)
    dof = np.array(dofs)
    mask = (z > -1e-6) & (tag == 0) & (dof <= 3)
    sel = az[mask]
    if sel.size == 0:
        raise RuntimeError('%s: zero free-surface non-fault non-PML nodes matched -- '
                            'mask is wrong or the dump is empty' % dump_path)
    return float(np.mean(np.abs(sel))), float(np.mean(sel)), int(sel.size)


def main():
    status = subprocess.run(['git', '-C', REPO_ROOT, 'status', '--porcelain',
                              '--', 'src/fortran/meshgen.f90'],
                             capture_output=True, text=True).stdout.strip()
    if status:
        raise SystemExit('FAIL: src/fortran/meshgen.f90 has uncommitted changes before this '
                          'script starts -- refusing to touch it:\n%s' % status)

    results = {}
    try:
        for label, present in (('before', True), ('after', False)):
            set_offset_present(present)
            rebuild()
            for case in CASES:
                case_dir = build_case(case)
                dump_path = run_dump(case_dir)
                mean_abs, mean_signed, n = stats(dump_path)
                results.setdefault(case, {})[label] = (mean_abs, mean_signed, n)
                print('%-6s %-14s mean|a_z|=%.4f  signed mean a_z=%.4f  n=%d'
                      % (label, case, mean_abs, mean_signed, n))
    finally:
        # Always restore the committed source and rebuild, whatever happened above.
        subprocess.run(['git', '-C', REPO_ROOT, 'checkout', '--', 'src/fortran/meshgen.f90'],
                        check=True)
        rebuild()
        status = subprocess.run(['git', '-C', REPO_ROOT, 'status', '--porcelain',
                                  '--', 'src/fortran/meshgen.f90'],
                                 capture_output=True, text=True).stdout.strip()
        if status:
            raise SystemExit('FAIL: src/fortran/meshgen.f90 did not restore cleanly: %r' % status)
        print('\nsrc/fortran/meshgen.f90 restored to the committed tree and bin/eqdyna rebuilt '
              'from it (verified clean via git status).')

    print('\n%-14s %14s %14s %14s %14s' % (
        'case', 'mean|az| before', 'mean|az| after', 'signed before', 'signed after'))
    for case in CASES:
        b = results[case]['before']
        a = results[case]['after']
        print('%-14s %14.4f %14.4f %14.4f %14.4f' % (case, b[0], a[0], b[1], a[1]))


if __name__ == '__main__':
    main()

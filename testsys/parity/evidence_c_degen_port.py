#! /usr/bin/env python3
"""
Regenerates, FROM SCRATCH, the C_degen wedge-degeneration MESH-port parity
check between a freshly-built Fortran binary and the Python port
(src/python/eqdyna/meshgen.py's build_elements/_build_elements_scalar,
porting src/fortran/library_degeneration.f90's wedge()/reorder() +
meshgen.f90:91-101's type-13 retag + meshgen.f90:763-767's checkIsOnFault
C_degen>3 branch), on test.tpv36 (dx=500, C_degen=par.dip=15).
"""
import datetime
import json
import os
import shutil
import socket
import subprocess
import sys
import tempfile
import time

import numpy as np

TESTSYS = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(TESTSYS))
PROBE_SRC = os.path.join(TESTSYS, 'probe_wedge_mesh.f90')

sys.path.insert(0, REPO_ROOT)
sys.path.insert(0, os.path.join(REPO_ROOT, 'testsys', 'e2e'))
sys.path.insert(0, os.path.join(REPO_ROOT, 'src', 'python'))
import run_e2e
from eqdyna import eqdyna3d, meshgen, readInputFiles

CASE_NAME = 'test.tpv36'


def provenance():
    def first_line(cmd):
        r = subprocess.run(cmd, capture_output=True, text=True)
        out = (r.stdout or r.stderr).splitlines()
        return out[0] if out else '(unavailable)'

    sha = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=REPO_ROOT,
                          capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(['git', 'status', '--porcelain'], cwd=REPO_ROOT,
                                 capture_output=True, text=True).stdout.strip())
    load1, load5, load15 = os.getloadavg()
    return dict(
        date_utc=datetime.datetime.utcnow().isoformat() + 'Z',
        host=socket.gethostname(), sha=sha, dirty=dirty,
        loadavg='%.2f %.2f %.2f' % (load1, load5, load15), ncpu=os.cpu_count(),
        gfortran=first_line(['gfortran', '--version']),
        mpif90=first_line(['mpif90', '--version']),
    )


def print_provenance(prov):
    print('')
    print('==== provenance ====')
    for k in ('date_utc', 'host', 'sha', 'dirty', 'loadavg', 'ncpu', 'gfortran', 'mpif90'):
        print('  %-9s: %s' % (k, prov[k]))


def build_fresh_fortran_and_probe():
    """Clean rebuild of src/fortran (MACHINE=ubuntu), then compile+link the
    probe program against the SAME freshly-built .o object files. No
    existing src/fortran/*.f90 file is modified; the probe is a brand-new
    program linked against the existing compiled units (see
    probe_wedge_mesh.f90's header for why it duplicates allocInit's body
    instead of linking eqdyna3d.o directly)."""
    src = os.path.join(REPO_ROOT, 'src', 'fortran')
    env = dict(os.environ, MACHINE='ubuntu')
    subprocess.run(['bash', '-c', 'rm -f *.o eqdyna'], cwd=src)
    t0 = time.time()
    r = subprocess.run(['make', 'MACHINE=ubuntu'], cwd=src, env=env,
                        capture_output=True, text=True)
    build_elapsed = time.time() - t0
    if r.returncode != 0:
        print(r.stdout[-4000:])
        print(r.stderr[-4000:])
        raise SystemExit('FAIL: make MACHINE=ubuntu exited %d' % r.returncode)
    print('built src/eqdyna fresh in %.1fs' % build_elapsed)

    probe_o = os.path.join(src, '_probe_wedge_mesh.o')
    probe_bin = os.path.join(src, '_probe_wedge_mesh')
    r = subprocess.run(['mpif90', '-c', '-fopenmp', '-ffree-line-length-none', '-O3',
                         '-I/usr/include', PROBE_SRC, '-o', probe_o],
                        cwd=src, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-4000:])
        print(r.stderr[-4000:])
        raise SystemExit('FAIL: compiling probe_wedge_mesh.f90 exited %d' % r.returncode)

    objs = sorted(f for f in os.listdir(src)
                  if f.endswith('.o') and f != 'eqdyna3d.o'
                  and f != os.path.basename(probe_o))
    r = subprocess.run(['mpif90', '-fopenmp', '-ffree-line-length-none', '-O3', probe_o]
                        + objs + ['-o', probe_bin, '-L/usr/lib/x86_64-linux-gnu',
                                  '-lnetcdf', '-lnetcdff'],
                        cwd=src, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-4000:])
        print(r.stderr[-4000:])
        raise SystemExit('FAIL: linking probe_wedge_mesh exited %d' % r.returncode)
    return os.path.join(src, 'eqdyna'), probe_bin


def run_probe(probe_bin, case_dir):
    dst = os.path.join(case_dir, os.path.basename(probe_bin))
    shutil.copy(probe_bin, dst)
    out_txt = os.path.join(case_dir, 'probe_wedge_mesh_out.txt')
    if os.path.exists(out_txt):
        os.remove(out_txt)
    t0 = time.time()
    r = subprocess.run(['mpirun', '-np', '1', dst], cwd=case_dir,
                        capture_output=True, text=True)
    elapsed = time.time() - t0
    if r.returncode != 0:
        print(r.stdout[-4000:])
        print(r.stderr[-4000:])
        raise SystemExit('FAIL: probe_wedge_mesh exited %d on %s' % (r.returncode, case_dir))
    if not os.path.isfile(out_txt):
        raise SystemExit('FAIL: %s was not written by the probe' % out_txt)
    print('Fortran probe (mesh-gen only) on %s: %.1fs' % (CASE_NAME, elapsed))
    return out_txt


def parse_probe_output(path):
    scalars = {}
    rows = {}
    in_elems = False
    with open(path) as f:
        for ln in f:
            s = ln.strip()
            if s == 'BEGIN_ELEMENTS':
                in_elems = True
                continue
            if s == 'END_ELEMENTS':
                in_elems = False
                continue
            if in_elems:
                parts = s.split()
                eid = int(parts[0])
                rows[eid] = (int(parts[1]), [int(x) for x in parts[2:10]],
                             [float(x) for x in parts[10:15]])
            else:
                parts = s.split()
                if len(parts) == 2:
                    scalars[parts[0]] = int(parts[1])
    return scalars, rows


def run_python_mesh(case_dir):
    params, g = readInputFiles.build_params(case_dir)
    material = readInputFiles.read_bmaterial(
        os.path.join(case_dir, 'bMaterial.txt'), g['nmat'], g['n2mat'])
    xline, yline, zline, pmlb, _ = meshgen.build_grid_lines(params)
    meshCoor, nftnd, nsmp = meshgen.build_node_coordinates(xline, yline, zline, params)

    t0 = time.time()
    conn, elem_type, mat, depth = meshgen.build_elements(
        xline, yline, zline, params, pmlb, nsmp, material, meshCoor)
    vec_elapsed = time.time() - t0

    t0 = time.time()
    conn_s, elem_type_s, mat_s, depth_s = meshgen._build_elements_scalar(
        xline, yline, zline, params, pmlb, nsmp, material, meshCoor)
    scalar_elapsed = time.time() - t0

    scalar_vs_vector_ok = (np.array_equal(conn, conn_s) and
                            np.array_equal(elem_type, elem_type_s) and
                            np.array_equal(mat, mat_s) and
                            np.array_equal(depth, depth_s))

    return dict(meshCoor=meshCoor, nftnd=nftnd, conn=conn, elem_type=elem_type, mat=mat,
                vec_elapsed=vec_elapsed, scalar_elapsed=scalar_elapsed,
                scalar_vs_vector_ok=scalar_vs_vector_ok)


def compare_mesh(fort_scalars, fort_rows, py):
    """Item (a)'s exact-match gate. Returns (ok, lines)."""
    lines = []
    ok = True
    checks = [
        ('totalNumOfElements', fort_scalars.get('totalNumOfElements'), py['conn'].shape[0]),
        ('totalNumOfNodes', fort_scalars.get('totalNumOfNodes'), py['meshCoor'].shape[0] - 1),
        ('nftnd1', fort_scalars.get('nftnd1'), py['nftnd']),
    ]
    for t in (1, 2, 11, 12, 13):
        checks.append(('count_type%d' % t, fort_scalars.get('count_type%d' % t),
                        int(np.sum(py['elem_type'] == t))))
    for name, fval, pval in checks:
        this_ok = (fval == pval)
        ok = ok and this_ok
        lines.append('%-20s fortran=%-10s python=%-10s %s'
                      % (name, fval, pval, 'OK' if this_ok else 'MISMATCH'))

    mismatches = []
    for eid, (etype_f, nodes_f, mat_f) in fort_rows.items():
        i = eid - 1
        etype_p = int(py['elem_type'][i])
        nodes_p = py['conn'][i].tolist()
        mat_p = py['mat'][i].tolist()
        row_ok = (etype_p == etype_f and nodes_p == nodes_f and
                  max(abs(a - b) for a, b in zip(mat_p, mat_f)) <= 1e-6)
        if not row_ok:
            mismatches.append((eid, etype_f, nodes_f, mat_f, etype_p, nodes_p, mat_p))
    ok = ok and not mismatches
    lines.append('per-element (type 11/12/13) rows checked=%d mismatches=%d'
                  % (len(fort_rows), len(mismatches)))
    for m in mismatches[:5]:
        lines.append('  MISMATCH elem=%d fortran=(type=%d,nodes=%r,mat=%r) '
                      'python=(type=%d,nodes=%r,mat=%r)' % m)
    lines.append('scalar-vs-vectorized (Python internal oracle check) on the FULL '
                 'real mesh: %s' % ('OK' if py['scalar_vs_vector_ok'] else 'MISMATCH'))
    ok = ok and py['scalar_vs_vector_ok']
    return ok, lines


def attempt_dynamic(case_dir):
    """Item (b): attempt a handful of numpy-backend dynamics steps. Expected
    to raise (see module docstring) -- caught and reported, not papered
    over. Returns (attempted_ok, message)."""
    try:
        eqdyna3d.run_case(case_dir, nsteps=5, verbose=False, backend='numpy')
        return True, ('UNEXPECTED: dynamic run succeeded (5 numpy steps) -- the '
                       'FEM-kernel wedge-degeneration guard did not fire; re-check '
                       'whether this mesh actually contains type-11/12 elements.')
    except NotImplementedError as exc:
        return False, 'raised NotImplementedError as expected: %s' % exc


def main():
    os.makedirs(os.path.join(TESTSYS, 'evidence_output'), exist_ok=True)
    prov = provenance()
    print_provenance(prov)

    tmp = tempfile.mkdtemp(prefix='evidence_c_degen_')
    try:
        case_dir = os.path.join(tmp, CASE_NAME)
        run_e2e.make_serial_case(CASE_NAME, case_dir, run_e2e.base_env())

        fortran_bin, probe_bin = build_fresh_fortran_and_probe()
        probe_out = run_probe(probe_bin, case_dir)
        fort_scalars, fort_rows = parse_probe_output(probe_out)

        py = run_python_mesh(case_dir)
        print('python vectorized build_elements: %.2fs; scalar oracle: %.2fs'
              % (py['vec_elapsed'], py['scalar_elapsed']))

        ok_a, lines_a = compare_mesh(fort_scalars, fort_rows, py)
        print('')
        print('==== item (a): MESH parity (test.tpv36, C_degen=15) ====')
        for ln in lines_a:
            print('  ' + ln)
        print('  ITEM (a) VERDICT: %s' % ('EXACT MATCH' if ok_a else 'MISMATCH'))

        dyn_attempted_ok, dyn_msg = attempt_dynamic(case_dir)
        print('')
        print('==== item (b): dynamic-step evidence (report-only, NOT gated) ====')
        print('  numpy backend, nsteps=5: %s' % dyn_msg)
        print('  jax backend: NOT run separately -- build_solver_state\'s wedge-'
              'element guard runs before any backend dispatch, so jax would hit '
              'the identical NotImplementedError for the identical reason.')

        out = dict(generated_utc=datetime.datetime.utcnow().isoformat() + 'Z',
                   provenance=prov, case=CASE_NAME,
                   item_a_ok=ok_a, item_a_lines=lines_a,
                   item_b_dynamic_attempted_ok=dyn_attempted_ok,
                   item_b_dynamic_message=dyn_msg)
        ts = datetime.datetime.utcnow().strftime('%Y%m%dT%H%M%SZ')
        out_path = os.path.join(TESTSYS, 'evidence_output', 'evidence_c_degen_%s.json' % ts)
        with open(out_path, 'w') as f:
            json.dump(out, f, indent=2, default=str)
        print('')
        print('Wrote %s' % out_path)
        print('')
        print('This script is REPORT-ONLY and always exits 0: item (a) mismatch '
              '(if any) is the thing to go tell a human about, not something this '
              'script itself fails a CI gate on -- it is not wired into '
              'testsys/run.py and test.tpv36/tpv37 are NOT added to '
              'testNameList.py/testsys/matrix.py by this script.')
        if not ok_a:
            print('')
            print('*** ITEM (a) MISMATCH -- see lines above. Do not treat this '
                  'run as a clean port. ***')
        return 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == '__main__':
    sys.exit(main())

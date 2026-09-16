#! /usr/bin/env python3
"""
Strong-scaling measurement (report-only; never a red/green gate).

Fortran MPI: test.tpv8 at 1/2/4/8/16/32 ranks, ranks pinned to distinct
cores (taskset), decomposition matched to rank count. Python (JAX-CPU):
the standalone solver at 1/2/4/8 XLA intra-op threads, pinned likewise.
One configuration at a time — scaling numbers are meaningless under
self-contention. Prints a provenance-stamped table (git SHA, host,
loadavg, compiler) per PROJECT_RULES.md rule 6.

Usage: python3 testsys/perf/run_scaling.py [--case test.tpv8] [--out table.txt]
Env: EQDYNAROOT, MACHINE, PATH with bin: and scripts:.
"""
import os, re, shutil, subprocess, sys, time, tempfile, json

ROOT = os.environ.get('EQDYNAROOT') or os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))))
CASE = 'test.tpv8'
DECOMP = {1: (1, 1, 1), 2: (2, 1, 1), 4: (2, 2, 1), 8: (2, 2, 2),
          16: (4, 2, 2), 32: (4, 4, 2)}
FORTRAN_RANKS = [1, 2, 4, 8, 16, 32]
PY_THREADS = [1, 2, 4, 8]


def sh(cmd, **kw):
    r = subprocess.run(cmd, shell=True, text=True, capture_output=True, **kw)
    if r.returncode != 0:
        raise RuntimeError(f'FAIL ({r.returncode}): {cmd}\n{r.stderr[-800:]}')
    return r


def make_case(dst, nx, ny, nz):
    if os.path.exists(dst):
        shutil.rmtree(dst)
    sh(f'create.newcase {dst} {CASE}')
    p = os.path.join(dst, 'user_defined_params.py')
    s = open(p).read()
    s = re.sub(r'par\.nx\s*,\s*par\.ny\s*,\s*par\.nz\s*=.*',
               f'par.nx, par.ny, par.nz = {nx}, {ny}, {nz}', s)
    s = re.sub(r'^par\.nx = .*', f'par.nx = {nx}', s, flags=re.M)
    s = re.sub(r'^par\.ny = .*', f'par.ny = {ny}', s, flags=re.M)
    s = re.sub(r'^par\.nz = .*', f'par.nz = {nz}', s, flags=re.M)
    open(p, 'w').write(s)
    sh('./case.setup', cwd=dst)


def run_fortran(work, n):
    nx, ny, nz = DECOMP[n]
    d = os.path.join(work, f'f{n}')
    make_case(d, nx, ny, nz)
    cores = ','.join(str(c) for c in range(n))
    t0 = time.time()
    sh(f'taskset -c {cores} mpirun -np {n} {ROOT}/bin/eqdyna', cwd=d)
    return time.time() - t0


def run_python(work, threads):
    d = os.path.join(work, f'p{threads}')
    make_case(d, 1, 1, 1)
    env = os.environ.copy()
    env['PYTHONPATH'] = os.path.join(ROOT, 'src', 'python')
    env['XLA_FLAGS'] = (f'--xla_cpu_multi_thread_eigen={"true" if threads > 1 else "false"} '
                        f'intra_op_parallelism_threads={threads}')
    env['OMP_NUM_THREADS'] = str(threads)
    cores = ','.join(str(c) for c in range(threads))
    t0 = time.time()
    r = subprocess.run(f'taskset -c {cores} python3 -m eqdyna .',
                       shell=True, cwd=d, env=env, text=True, capture_output=True)
    if r.returncode != 0:
        raise RuntimeError(f'python run failed:\n{r.stderr[-800:]}')
    return time.time() - t0


def main():
    global CASE
    if '--case' in sys.argv:
        CASE = sys.argv[sys.argv.index('--case') + 1]
    sha = sh(f'git -C {ROOT} rev-parse --short HEAD').stdout.strip()
    host = os.uname().nodename
    load = os.getloadavg()
    work = tempfile.mkdtemp(prefix='scaling.', dir=os.environ.get('TMPDIR', '/tmp'))
    rows = []
    tf1 = None
    for n in FORTRAN_RANKS:
        t = run_fortran(work, n)
        tf1 = tf1 or t
        rows.append(('fortran', n, t, tf1 / t, tf1 / t / n))
        print(f'fortran np={n:<3} {t:8.1f}s  speedup {tf1/t:5.2f}x  efficiency {tf1/t/n:5.1%}', flush=True)
    tp1 = None
    for n in PY_THREADS:
        t = run_python(work, n)
        tp1 = tp1 or t
        rows.append(('python-jax', n, t, tp1 / t, tp1 / t / n))
        print(f'python  th={n:<3} {t:8.1f}s  speedup {tp1/t:5.2f}x  efficiency {tp1/t/n:5.1%}', flush=True)
    meta = dict(case=CASE, sha=sha, host=host, loadavg=load,
                date=time.strftime('%Y-%m-%d %H:%M'), rows=rows)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'scaling_last.json')
    json.dump(meta, open(out, 'w'), indent=1)
    print(f'provenance: {CASE} @ {sha} on {host}, loadavg {load}; saved {out}')
    shutil.rmtree(work, ignore_errors=True)


if __name__ == '__main__':
    main()

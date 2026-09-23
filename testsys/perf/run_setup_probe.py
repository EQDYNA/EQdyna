#! /usr/bin/env python3
"""STAGE 0 of the setup-rewrite question (pathway item 64): how expensive IS
the per-rank setup, and does it get worse when 32 ranks do it at once?

WHAT THIS MEASURES, AND WHY IT IS NOT A SOLVER RUN. `eqdyna3d.build_solver_state`
is the whole serial mesh build + input read that EVERY python-MPI rank runs in
full before the solve begins (`MPI4NodalQuant.py:35-41`: the port "keeps the
PROVEN serial mesh build on every rank and then RESTRICTS it"). Item 64 suspects
that 32 concurrent copies of this, not the stepping, is what makes the 32-rank
WALL per-step read 121.73 ms against a differenced 45.76 ms. This probe runs
build_solver_state ALONE -- no time loop, no jax, no MPI -- at 1 process and at
N concurrent processes, and reports per-process wall, CPU-seconds, effective
cores and peak RSS.

TWO PREMISES THIS PROBE EXISTS TO SETTLE, both of which were WRONG in the
briefing that commissioned the setup rewrite (corrected in item 64 on
2026-09-22, from source):
  * the Fortran ALREADY builds rank-locally (`meshgen.f90:65-67` loops
    `do ix = 1, nx` with nx rank-local), so the global mesh build is
    PYTHON-ONLY. Nothing here measures the Fortran.
  * the "~100 s setup" figure has NEVER been measured by anything:
    `writeCompTime` is declared 0 at `globalvar.f90:109` and set to 1 nowhere,
    so `compTimeInSeconds(1)` has never printed. Until this probe runs there is
    no setup-time number for this code at all.

MEASURED, NOT PROBED. Effective cores per process is (utime+stime)/wall from
the child's OWN getrusage after the build -- not a /proc/stat reading taken
before launch, which on this box has read busy 0.00 for cpus that then
delivered a quarter core each (cpus 0, 1 and 16-19 on record).

PARTIAL-LOSS SAFE (rule 20a). Every child writes its own JSON the moment it
finishes and the parent rewrites the snapshot after every repetition, so a
killed run keeps everything already measured.
"""
import argparse
import json
import os
import subprocess
import sys
import time

TESTSYS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.environ.get('EQDYNAROOT') or os.path.dirname(TESTSYS)
PKG = os.path.join(ROOT, 'src', 'python')
sys.path.insert(0, os.path.join(TESTSYS, 'perf'))

CHILD = r'''
import json, os, resource, sys, time
sys.path.insert(0, %(pkg)r)
cpu = %(cpu)d
if cpu >= 0:
    os.sched_setaffinity(0, {cpu})
t0 = time.perf_counter()
from eqdyna import eqdyna3d
t_import = time.perf_counter() - t0
r0 = resource.getrusage(resource.RUSAGE_SELF)
t1 = time.perf_counter()
S, mesh = eqdyna3d.build_solver_state(%(case)r)
t_build = time.perf_counter() - t1
r1 = resource.getrusage(resource.RUSAGE_SELF)
cpu_s = (r1.ru_utime - r0.ru_utime) + (r1.ru_stime - r0.ru_stime)
out = dict(idx=%(idx)d, cpu=cpu, import_s=t_import, build_s=t_build,
           build_cpu_s=cpu_s, eff_cores=cpu_s / t_build if t_build else None,
           maxrss_mb=r1.ru_maxrss / 1024.0,
           nodes=int(mesh['meshCoor'].shape[-1]) if hasattr(mesh.get('meshCoor'), 'shape') else None)
open(%(dst)r, 'w').write(json.dumps(out))
print(json.dumps(out), flush=True)
'''


def cpu_pool(exclude):
    """CPUs this process may use, minus the ones MEASURED to misbehave."""
    return [c for c in sorted(os.sched_getaffinity(0)) if c not in exclude]


def one_wave(case_dir, n, cpus, outdir, tag, backend_env):
    """n concurrent build_solver_state processes, one pinned per cpu."""
    procs, dsts = [], []
    env = dict(os.environ)
    env.update(backend_env)
    t0 = time.time()
    for i in range(n):
        dst = os.path.join(outdir, '%s_%d.json' % (tag, i))
        dsts.append(dst)
        src = CHILD % dict(pkg=PKG, cpu=cpus[i % len(cpus)], case=case_dir,
                           idx=i, dst=dst)
        procs.append(subprocess.Popen([sys.executable, '-c', src], env=env,
                                      stdout=subprocess.DEVNULL,
                                      stderr=subprocess.PIPE))
    errs = []
    for p in procs:
        _, e = p.communicate()
        if p.returncode != 0:
            errs.append(e.decode()[-400:])
    wall = time.time() - t0
    rows = []
    for d in dsts:
        if os.path.isfile(d):
            rows.append(json.load(open(d)))
    return dict(n=n, wave_wall_s=wall, rows=rows, failures=errs,
                completed=len(rows))


def summarise(w):
    if not w['rows']:
        return '  n=%-3d ALL %d PROCESSES FAILED (wave wall %.1f s)' % (
            w['n'], w['n'], w['wave_wall_s'])
    b = [r['build_s'] for r in w['rows']]
    e = [r['eff_cores'] for r in w['rows'] if r['eff_cores']]
    m = [r['maxrss_mb'] for r in w['rows']]
    i = [r['import_s'] for r in w['rows']]
    return ('  n=%-3d build_s min/mean/max %.2f/%.2f/%.2f   import_s mean %.2f   '
            'eff cores min/mean %.2f/%.2f   maxrss MB mean %.0f (sum %.0f)   '
            'wave wall %.1f s   %d/%d completed'
            % (w['n'], min(b), sum(b) / len(b), max(b), sum(i) / len(i),
               min(e) if e else -1, sum(e) / len(e) if e else -1,
               sum(m) / len(m), sum(m), w['wave_wall_s'], w['completed'], w['n']))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', default='test.tpv104')
    ap.add_argument('--procs', default='1,32')
    ap.add_argument('--reps', type=int, default=3)
    ap.add_argument('--exclude-cpus', default='0,1,16,17,18,19')
    ap.add_argument('--case-dir', default='',
                    help='reuse an already-built serial case dir instead of '
                         'building one (build_solver_state only READS it)')
    ap.add_argument('--out', default='')
    a = ap.parse_args()
    ns = [int(x) for x in a.procs.split(',') if x]
    excl = {int(x) for x in a.exclude_cpus.split(',') if x.strip()}
    cpus = cpu_pool(excl)
    if len(cpus) < max(ns):
        raise SystemExit('FAIL: %d usable cpus after exclusions, need %d'
                         % (len(cpus), max(ns)))

    case_dir = a.case_dir
    if not case_dir:
        import run_scaling as rs
        rs.CASE = a.case
        case_dir = rs.build_py_case(a.case)
    if not os.path.isdir(case_dir):
        raise SystemExit('FAIL: no case dir %s' % case_dir)

    out = a.out or os.path.join(
        ROOT, 'docs', 'perf_snapshots',
        'setup_probe_%s.json' % time.strftime('%Y-%m-%d_%H%M%S'))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    outdir = os.path.join(os.path.dirname(out),
                          'setup_probe_children_%s' % time.strftime('%H%M%S'))
    os.makedirs(outdir, exist_ok=True)

    sha = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=ROOT,
                         capture_output=True, text=True).stdout.strip()
    print('## setup probe  case %s  sha %s  cpus %d usable (excluded %s)  '
          'loadavg %.2f' % (a.case, sha, len(cpus), sorted(excl),
                            os.getloadavg()[0]), flush=True)
    print('   case dir %s' % case_dir, flush=True)

    waves = []

    def snap():
        json.dump(dict(case=a.case, case_dir=case_dir, sha=sha,
                       host=os.uname().nodename,
                       date=time.strftime('%Y-%m-%d %H:%M'),
                       excluded_cpus=sorted(excl), usable_cpus=cpus,
                       reps=a.reps, waves=waves,
                       what='eqdyna3d.build_solver_state only -- no time loop, '
                            'no MPI, no jax'),
                  open(out, 'w'), indent=1)

    for rep in range(1, a.reps + 1):
        for n in ns:
            print('-- rep %d  n=%d  loadavg %.2f' % (rep, n, os.getloadavg()[0]),
                  flush=True)
            w = one_wave(case_dir, n, cpus, outdir, 'r%d_n%d' % (rep, n),
                         {'EQDYNA_BACKEND': 'numpy'})
            w['rep'] = rep
            w['loadavg'] = os.getloadavg()
            waves.append(w)
            print(summarise(w), flush=True)
            for f in w['failures'][:2]:
                print('  CHILD STDERR: %s' % f.replace('\n', ' | '), flush=True)
            snap()

    print('\nsnapshot %s' % out, flush=True)
    for n in ns:
        got = [w for w in waves if w['n'] == n and w['rows']]
        if not got:
            continue
        per = [r['build_s'] for w in got for r in w['rows']]
        print('SUMMARY n=%-3d build_s mean %.2f over %d process-samples in %d '
              'wave(s)' % (n, sum(per) / len(per), len(per), len(got)),
              flush=True)


if __name__ == '__main__':
    main()

#! /usr/bin/env python3
"""
Per-step cost vs elements-per-rank (report-only; never a red/green gate).

Tests the claim "EQdyna's speed is a linear function of cells per rank" by
crossing resolution with rank count so the two are decorrelated: the same
elements-per-rank value is reached from several (dx, ranks) combinations.

Design decisions and why:

* npy = 1 in every decomposition. EQdyna's fault is the x-z plane at y=0, and
  a y partition boundary landing on it is a known defect class
  (pathway_forward.md, "fault-plane-on-MPI-boundary halving"). npy=1 removes
  the condition entirely rather than relying on the fix.
* Fixed, short step count per configuration (default 200). Per-step cost is
  the measurand; physics is not. `nstep = idnint(term/dt)`
  (readInputFiles.f90:128) and `dt = 0.5*dx/vp`, so `term` is set per dx.
* Seconds per timestep is taken from the run's OWN clocks, not from wall time
  alone. `library_output.f90:output_timeanalysis` writes `compTime<me>` for
  every rank with MPI_WTIME totals:
      1 setup+mesh, 2 mass assembly, 3 velDispUpdate, 4 assembleGlobalKU,
      5 calcHourglassResist, 6 faulting, 8 output, 9 whole program
  so the time-stepping loop is 9 - 1 - 8 (comp(2) is NOT subtracted: it is
  corrupted -- see the comment at the `loop =` line), and per-step is that
  over nstep. `kernel` (comp 3..6, accumulated only inside driver.f90's loop)
  and `halo` (MPICommTime) are recorded separately, which is what separates
  element-proportional work from rank-count-dependent overhead. Wall time is
  recorded alongside as an independent cross-check, as is the last
  `TimeElapsed (s)` value the run printed -- it equals `term` when nstep steps
  ran. (The NUMBER of those lines is not the step count:
  faulting.f90:showSourceDynamics fires once per fault node pair matching the
  hypocentre, which is 1, 2 or 4 pairs depending on dx.)
  This needs `writeCompTime = 1` (globalvar.f90:109), which is off in the
  default build -- point --bin at a build that has it on.
* Elements per rank comes from the same `compTime<me>` files (each rank
  writes its own `totalNumOfElements`), cross-checked against the standalone
  replica in elem_per_rank.py. Neither the "GB memory is expected" log line
  nor case.setup's estimate is a per-rank element count -- see elem_per_rank.py.

Machine courtesy: one configuration at a time, and each launch waits for the
box to be free of other eqdyna processes and for the load average to drop
(--max-others / --load-ceiling). Every row records the load it started under,
whether another eqdyna job appeared, and the max/min spread of per-rank loop
time -- a row measured while a core was shared shows up in that spread.

Usage:
  python3 testsys/perf/run_elem_scaling.py --bin <eqdyna> --work <dir> [--out f.json]
  python3 testsys/perf/run_elem_scaling.py ... --grid 500:16,250:5x1x1,125:5x1x5
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import elem_per_rank  # noqa: E402

CASE = 'test.tpv8'
# npy is 1 everywhere: the fault plane is never on a partition boundary.
DECOMP = {1: (1, 1, 1), 4: (2, 1, 2), 16: (4, 1, 4), 48: (8, 1, 6)}
MAX_RANKS = 48


def sh(cmd, **kw):
    r = subprocess.run(cmd, shell=True, text=True, capture_output=True, **kw)
    if r.returncode != 0:
        raise RuntimeError(f'FAIL ({r.returncode}): {cmd}\n{r.stdout[-400:]}\n{r.stderr[-800:]}')
    return r


def others_running():
    r = subprocess.run('pgrep -c eqdyna', shell=True, text=True, capture_output=True)
    return int(r.stdout.strip() or 0)


def wait_for_idle(need_ranks, poll=60, load_ceiling=12.0, skip=False,
                  max_others=0):
    """Block until the box is quiet enough to measure on.

    `max_others` is normally 0. It exists because this box is shared: other
    agents run short 4-rank cases here, and a strict zero can stay closed
    indefinitely. Raising it trades a guaranteed-quiet box for a flagged one --
    every row records `contended` and `loop_spread`, and a row whose ranks
    actually had to share cores shows a large `loop_spread`.
    """
    if skip:
        return os.getloadavg()[0]
    clear = 0
    while True:
        n, load = others_running(), os.getloadavg()[0]
        if n <= max_others and load < load_ceiling:
            clear += 1
            # two consecutive clear polls: the owner's queue can start the next
            # job seconds after the previous one exits, and a single clear
            # sample would launch us straight into it.
            if clear >= 2:
                return load
            print(f'  [wait] clear ({load:.1f}); confirming in {poll}s', flush=True)
        else:
            clear = 0
            print(f'  [wait] {n} eqdyna proc(s), load {load:.1f}; need <= '
                  f'{max_others} procs and load < {load_ceiling}', flush=True)
        time.sleep(poll)


def make_case(dst, dx, npx, npy, npz, nstep):
    if os.path.exists(dst):
        shutil.rmtree(dst)
    env = dict(os.environ, EQDYNAROOT=ROOT)
    subprocess.run(f'{ROOT}/scripts/create.newcase {dst} {CASE}', shell=True,
                   check=True, env=env, capture_output=True, text=True)
    p = os.path.join(dst, 'user_defined_params.py')
    s = open(p).read()
    s = re.sub(r'^par\.dx = .*$', f'par.dx = {float(dx)}', s, flags=re.M)
    s = re.sub(r'^par\.nx = .*$', f'par.nx = {npx}', s, flags=re.M)
    s = re.sub(r'^par\.ny = .*$', f'par.ny = {npy}', s, flags=re.M)
    s = re.sub(r'^par\.nz = .*$', f'par.nz = {npz}', s, flags=re.M)
    # dt is set from dx further down the file; pin term so nstep is exactly nstep.
    s = re.sub(r'^par\.term = .*$', f'par.term = {nstep}*0.5*par.dx/par.vp', s, flags=re.M)
    open(p, 'w').write(s)
    sh('./case.setup', cwd=dst, env=env)


def read_comptime(d, nranks):
    """Per-rank timing + element count from compTime<me> (library_output.f90:224)."""
    rows = []
    for me in range(nranks):
        f = os.path.join(d, f'compTime{me}')
        if not os.path.exists(f):
            raise RuntimeError(f'missing {f} -- is writeCompTime=1 in this build?')
        v = open(f).read().split()
        comp = [float(x) for x in v[:9]]
        rows.append(dict(me=me, comp=comp, mpicomm=float(v[9]),
                         nelem=int(v[10]), neq=int(v[11])))
    return rows


def parse_cfg(tok):
    """'500:16' -> (500, (4,1,4)); '250:5x1x1' -> (250, (5,1,1))."""
    dxs, decs = tok.split(':')
    if 'x' in decs:
        dec = tuple(int(v) for v in decs.split('x'))
    else:
        dec = DECOMP[int(decs)]
    return int(dxs), dec


def run_one(work, dx, dec, nstep, binary, extra_mpi='', skip_wait=False,
            load_ceiling=12.0, max_others=0):
    npx, npy, npz = dec
    nranks = npx * npy * npz
    d = os.path.join(work, f'dx{dx}_np{nranks}_{npx}x{npy}x{npz}')
    print(f'[setup ] dx={dx} ranks={nranks} decomp={npx}x{npy}x{npz} nstep={nstep}', flush=True)
    make_case(d, dx, npx, npy, npz, nstep)

    pred = elem_per_rank.counts(elem_per_rank.read_case(d), npx, npy, npz)
    print(f'         predicted elem/rank max {pred["elem_max"]:,} mean {pred["elem_mean"]:,.0f}',
          flush=True)

    load_before = wait_for_idle(nranks, load_ceiling=load_ceiling, skip=skip_wait,
                                max_others=max_others)
    env = dict(os.environ, OMP_NUM_THREADS='1', EQDYNAROOT=ROOT)
    # Hard 1:1 pin, ranks on cores 0..nranks-1, the same policy for every
    # configuration so NUMA locality is a known function of rank count only.
    # Verified on this box (Open MPI 4.1.1, checked with Cpus_allowed_list in
    # the children): `--cpu-set` alone confines all ranks to the set but does
    # NOT pin them 1:1; a `taskset` mask on mpirun is simply overridden (Open
    # MPI calls sched_setaffinity itself and lands on cores 0..n-1 anyway);
    # and `--cpu-set` with `--map-by core` is rejected outright ("Conflicting
    # directives for mapping policy"). Only this form gives one rank per core.
    cmd = (f'mpirun --bind-to core --map-by core --report-bindings '
           f'{extra_mpi} -np {nranks} {binary}')
    t0 = time.time()
    r = subprocess.run(cmd, shell=True, cwd=d, env=env, text=True, capture_output=True)
    wall = time.time() - t0
    contended = others_running()  # >0 means someone else's job started mid-run
    open(os.path.join(d, 'run.stdout'), 'w').write(r.stdout)
    open(os.path.join(d, 'run.stderr'), 'w').write(r.stderr)
    if r.returncode != 0:
        raise RuntimeError(f'run failed dx={dx} np={nranks}:\n{r.stdout[-600:]}\n{r.stderr[-800:]}')

    rows = read_comptime(d, nranks)
    # nstep as the binary actually saw it, from the run's own log.
    def logval(label):
        m = re.search(label + r'\s*:\s*([0-9.E+-]+)', r.stdout)
        return float(m.group(1)) if m else None
    log_cells_rank0 = logval('Cells per rank')
    log_cells_total = logval('Cells total')
    te = re.findall(r'TimeElapsed\s+\(s\)\s+([0-9.E+-]+)', r.stdout)
    # comp(9) - comp(1) - comp(8): whole program, less setup/mesh, less output.
    # comp(2) (mass assembly) is deliberately NOT subtracted -- it is corrupted:
    # assembleGlobalMass -> MPI4NodalQuant resets the SHARED global
    # `startTimeStamp` (assembleGlobalMass.f90:67), so eqdyna3d.f90:69 measures
    # only the tail of the call. The real mass-assembly cost (a one-off,
    # O(elements), roughly one timestep's worth) therefore stays inside `loop`,
    # biasing per-step cost by ~1/nstep and, being element-proportional, adding
    # no spurious intercept. `kernel` below is the uncontaminated inner-loop
    # number: comp(3..6) are accumulated ONLY from inside driver.f90's loop.
    loop = [x['comp'][8] - x['comp'][0] - x['comp'][7] for x in rows]
    kernel = [x['comp'][2] + x['comp'][3] + x['comp'][4] + x['comp'][5] for x in rows]
    p = elem_per_rank.read_case(d)
    nstep_actual = round(p['term'] / p['dt'])
    rec = dict(
        dx=dx, ranks=nranks, decomp=[npx, npy, npz], nstep=nstep_actual,
        dt=p['dt'], term=p['term'],
        elem_max=max(x['nelem'] for x in rows),
        elem_min=min(x['nelem'] for x in rows),
        elem_mean=sum(x['nelem'] for x in rows) / nranks,
        elem_rank0=rows[0]['nelem'],
        pred_elem_max=pred['elem_max'], pred_elem_mean=pred['elem_mean'],
        pred_per_rank_match=sorted(x['nelem'] for x in rows) == sorted(pred['per_rank']),
        global_elements=pred['global_elements'],
        wall=wall,
        t_total_max=max(x['comp'][8] for x in rows),
        t_setup_max=max(x['comp'][0] for x in rows),
        t_mass_max=max(x['comp'][1] for x in rows),
        t_output_max=max(x['comp'][7] for x in rows),
        t_loop_max=max(loop), t_loop_min=min(loop),
        t_loop_mean=sum(loop) / nranks,
        # loop_spread >> elem_spread means a rank lost time it should not have
        # (a core shared with someone else's job, or OS jitter) -- the check
        # that a row measured on a shared box is still usable.
        loop_spread=max(loop) / min(loop) - 1.0,
        elem_spread=max(q['nelem'] for q in rows) / min(q['nelem'] for q in rows) - 1.0,
        t_kernel_max=max(kernel),
        sec_per_step=max(loop) / nstep_actual,
        sec_per_step_kernel=max(kernel) / nstep_actual,
        sec_per_step_mpi=max(x['mpicomm'] for x in rows) / nstep_actual,
        sec_per_step_wall=wall / nstep_actual,
        t_veldisp=max(x['comp'][2] for x in rows),
        t_ku=max(x['comp'][3] for x in rows),
        t_hg=max(x['comp'][4] for x in rows),
        t_fault=max(x['comp'][5] for x in rows),
        t_mpicomm=max(x['mpicomm'] for x in rows),
        log_cells_rank0=log_cells_rank0, log_cells_total=log_cells_total,
        last_TimeElapsed=float(te[-1]) if te else None,
        n_TimeElapsed_lines=len(te),
        load_before=load_before, load_after=os.getloadavg()[0],
        contended=contended,
        percore=[x['nelem'] for x in rows],
    )
    print(f'[done  ] elem/rank {rec["elem_max"]:,} (max)  loop {rec["t_loop_max"]:.2f}s / '
          f'{nstep_actual} steps = {rec["sec_per_step"]*1e3:.2f} ms/step  '
          f'(wall {wall:.1f}s)  replica-match={rec["pred_per_rank_match"]}'
          + f'  spread {rec["loop_spread"]:+.1%} (elem {rec["elem_spread"]:+.1%})'
          + ('  *** OTHER EQDYNA RUNNING ***' if contended else ''), flush=True)
    # keep only the small text outputs; the case tree is regenerable
    for junk in ('on_fault_vars_input.nc',):
        f = os.path.join(d, junk)
        if os.path.exists(f):
            os.remove(f)
    return rec


def provenance():
    def out(c):
        return subprocess.run(c, shell=True, text=True, capture_output=True).stdout.strip()
    return dict(
        sha=out(f'git -C {ROOT} rev-parse --short HEAD'),
        dirty=bool(out(f'git -C {ROOT} status --porcelain -- src')),
        host=os.uname().nodename, kernel=os.uname().release,
        cpu=out("lscpu | grep 'Model name' | cut -d: -f2").strip(),
        cores=os.cpu_count(),
        fc=out('mpif90 --version | head -1'),
        mpi=out('mpirun --version | head -1'),
        netcdf=out('nc-config --version'),
        machine=os.environ.get('MACHINE', 'ubuntu'),
        fflags='-fopenmp -ffree-line-length-none -O3 (src/makefile, MACHINE=ubuntu)',
        date=time.strftime('%Y-%m-%d %H:%M:%S %Z'),
        case=CASE,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bin', required=True)
    ap.add_argument('--work', required=True)
    ap.add_argument('--out', default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                  'elem_scaling_last.json'))
    ap.add_argument('--nstep', type=int, default=200)
    ap.add_argument('--grid', default='500:1,500:4,500:16,500:48,'
                                      '250:1,250:4,250:16,250:48,'
                                      '125:4,125:16,125:48')
    ap.add_argument('--load-ceiling', type=float, default=12.0,
                    help='max 1-min load average to launch at. This box also '
                         'runs two long GNS/PyTorch trainings that hold the '
                         'idle load near 7, so a lower ceiling never opens.')
    ap.add_argument('--max-others', type=int, default=0,
                    help='how many other eqdyna processes may already be on the '
                         'box at launch (default 0). Raise it only to coexist '
                         'with small sibling jobs; every row still records '
                         'contended/loop_spread so a spoilt row is visible.')
    ap.add_argument('--smoke', action='store_true',
                    help='harness self-check only: skip the machine-courtesy wait '
                         '(use ONLY for a 1-rank, few-step validation run)')
    ap.add_argument('--nstep-override', default='',
                    help='dx:ranks:nstep,... for configs too expensive at --nstep')
    a = ap.parse_args()

    override = {}
    for tok in filter(None, a.nstep_override.split(',')):
        dx, nr, ns = tok.split(':')
        override[(int(dx), int(nr))] = int(ns)

    os.makedirs(a.work, exist_ok=True)
    recs = []
    if os.path.exists(a.out):
        recs = json.load(open(a.out)).get('rows', [])
    done = {(r['dx'], r['ranks']) for r in recs}
    prov = provenance()
    print(json.dumps(prov, indent=1), flush=True)

    for tok in a.grid.split(','):
        dx, dec = parse_cfg(tok)
        nranks = dec[0] * dec[1] * dec[2]
        if nranks > MAX_RANKS:
            raise SystemExit(f'refusing {nranks} ranks (cap {MAX_RANKS})')
        if dec[1] != 1:
            raise SystemExit(f'refusing npy={dec[1]}: npy must be 1 so no MPI '
                             'partition boundary can land on the fault plane')
        if (dx, nranks) in done:
            print(f'[skip  ] dx={dx} ranks={nranks} already in {a.out}', flush=True)
            continue
        rec = run_one(a.work, dx, dec, override.get((dx, nranks), a.nstep),
                      a.bin, skip_wait=a.smoke, load_ceiling=a.load_ceiling,
                      max_others=a.max_others)
        recs.append(rec)
        json.dump(dict(provenance=prov, rows=recs), open(a.out, 'w'), indent=1)
    print(f'wrote {a.out} ({len(recs)} rows)')


if __name__ == '__main__':
    main()

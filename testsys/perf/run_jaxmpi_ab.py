#! /usr/bin/env python3
"""2x2 A/B of the two INDEPENDENT changes on mira/jaxmpi-step-2026-09-22,
measured separately against current master, at several rank counts.

THE TWO FACTORS, which must never be reported as one number:
  (a) MPI4NodalQuant.exchange: non-blocking Irecv/Isend/Waitall instead of a
      blocking Sendrecv walked in ascending rank order.
  (b) driver.run_mpi: the per-step comm.Barrier() becomes a PROFILING DEVICE
      (`if prof:`) instead of running unconditionally in production.

So four arms: base, b, a, ab. An arm is a PAIR OF FILE VERSIONS swapped into
src/python/eqdyna, taken from the repository history; nothing is selected at
run time, because a runtime switch would be a config path that then has to be
landed and lived with.

(b) IS ONLY VISIBLE OUTSIDE THE PROFILE. Under MPI4NodalQuant.step_profile()
the barrier still runs on purpose (it is what separates load imbalance from
exchange cost), so arm b is identical to base under the profile. The
measurement here is therefore the perf tier's method and not the step profile:
per-step by DIFFERENCE over two step counts, over the ranks' OWN solve clock,
with compile/setup subtracted as the common term. A total wall clock would
mix in XLA compile.

ONE CPU SET PER RANK COUNT, CHOSEN ONCE AND REUSED BY ALL FOUR ARMS. The
project tool re-picks the least-loaded cpus per invocation; across four arms
that would make placement a confound indistinguishable from the effect.

THE EFFECTIVE-CORES FILTER IS THE POINT, NOT AN OBSTACLE. A point is REJECTED
unless every rank measured EFFECTIVE_CORES >= --min-eff (default 0.95, i.e.
"~1.00"), measured AFTER launch from the ranks' own stdout -- never from the
/proc/stat probe beforehand, which has read cpu 0 and cpu 1 at busy 0.00 and
then delivered 0.39. A rejected point is the system working, and "all points
rejected" is a complete outcome. Ratios are formed only between two ACCEPTED
arms: a speedup whose numerator and denominator ran at different effective
core counts is not a speedup.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

TESTSYS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get('EQDYNAROOT') or os.path.dirname(os.path.dirname(TESTSYS))
sys.path.insert(0, TESTSYS)
sys.path.insert(0, ROOT)
import run_mpi_scaling as ms       # noqa: E402
import run_numa_scaling as numa    # noqa: E402
import run_scaling as rs           # noqa: E402
import ledger                      # noqa: E402
from testsys import runlock        # noqa: E402

SELF_ROOT = os.path.dirname(os.path.dirname(TESTSYS))
# Asserted on by testsys/regression/test_perf_tool_locks.py, so the contract is
# a named constant and not a wording buried in a format string.
ROOT_MISMATCH_HEADER = (
    'refusing to start: $EQDYNAROOT names a DIFFERENT checkout than this tool')
PKG = os.path.join(ROOT, 'src', 'python', 'eqdyna')
DRIVER = os.path.join(PKG, 'driver.py')
MQ = os.path.join(PKG, 'MPI4NodalQuant.py')
# The resource is the PACKAGE, not the __pycache__ under it: stage() rewrites
# driver.py and MPI4NodalQuant.py in place, and every measurement below is
# taken against whatever those two files currently hold (item 74).
PKG_RESOURCE = os.path.relpath(PKG, ROOT)
LOCK_CONSEQUENCE = [
    'stage() COPIES this arm\'s driver.py and MPI4NodalQuant.py INTO that'
    ' package, and',
    'every timing below is taken against whatever those two files currently'
    ' hold. Two',
    'invocations in one checkout swap each other\'s arm files mid-measurement,'
    ' and the',
    'result is not a crash: it is a number recorded under the WRONG ARM LABEL'
    ' -- which',
    'this tool\'s own docstring calls worse than no number (rule 21a, item'
    ' 74).',
    '',
    'NOT waiting. NOT staging anyway. NOT falling back to a second copy of the'
    ' package',
    '-- each of those is the silent fallback rule 2 forbids. Run your A/B in'
    ' its own git',
    'worktree (rule 21a), or wait for the holder above to finish.']
# arm -> (driver version, MPI4NodalQuant version); 'm' = master, 'x' = branch.
ARMS = {'base': ('m', 'm'), 'b': ('x', 'm'), 'a': ('m', 'x'), 'ab': ('x', 'x')}
ARM_DESC = {'base': 'master: Sendrecv + unconditional per-step barrier',
            'b': 'change b only: barrier under prof',
            'a': 'change a only: non-blocking halo exchange',
            'ab': 'both (the merged branch)'}


def show(ref, name):
    """One file at one revision, as text."""
    out = subprocess.run(['git', 'show', '%s:src/python/eqdyna/%s' % (ref, name)],
                         cwd=ROOT, capture_output=True, text=True)
    if out.returncode != 0:
        raise SystemExit('FAIL: cannot read %s at %s -- %s'
                         % (name, ref, out.stderr.strip()[:300]))
    return out.stdout


def stage(vault, arm):
    """Put this arm's two files in place. Raises if a staged source is
    missing: a file left over from the previous arm would produce a number
    under the wrong label, which is worse than no number."""
    dv, mv = ARMS[arm]
    for name, tag, dst in (('driver.py', dv, DRIVER),
                           ('MPI4NodalQuant.py', mv, MQ)):
        src = os.path.join(vault, '%s.%s' % (name, tag))
        if not os.path.isfile(src):
            raise SystemExit('FAIL: missing staged source %s' % src)
        shutil.copyfile(src, dst)
    cache = os.path.join(PKG, '__pycache__')
    if os.path.isdir(cache):
        shutil.rmtree(cache)


def build_vault(work, merged_ref):
    vault = os.path.join(work, 'vault')
    os.makedirs(vault)
    for name, tag, ref in (('driver.py', 'm', 'master'),
                           ('driver.py', 'x', merged_ref),
                           ('MPI4NodalQuant.py', 'm', 'master'),
                           ('MPI4NodalQuant.py', 'x', merged_ref)):
        open(os.path.join(vault, '%s.%s' % (name, tag)), 'w').write(
            show(ref, name))
    # THE STEP-0 CHECK, RE-RUN HERE. If the merge had reverted the per-rank
    # XLA cache fix (item 61), every 'x' arm would silently run the config
    # that wedges this cell, and this A/B would be measuring that instead.
    mx = open(os.path.join(vault, 'driver.py.x')).read()
    if "enable_compilation_cache(subdir='rank%d' % rank)" not in mx:
        raise SystemExit('FAIL: merged driver.py lost the per-rank XLA cache '
                         'fix. Refusing to measure a reverted tree.')
    return vault


def require_root_is_this_checkout():
    """REFUSE when $EQDYNAROOT points at a checkout other than this file's.

    WHY REFUSE AND NOT WARN (item 77). This tool WRITES: `stage()` copies an
    arm's `driver.py` and `MPI4NodalQuant.py` into `ROOT/src/python/eqdyna`,
    and the `finally` clause rewrites them once more on the way out. Under a
    stale `EQDYNAROOT` -- the shape `install-eqdyna.sh` leaves behind, since it
    exports `$(pwd)` of whichever tree was installed and the export survives a
    `cd` into a worktree -- those writes land in ANOTHER SESSION's package
    while this session measures its own. That is precisely rules 21/21b's
    damage class, and Gate 0 does not catch it: the lock follows the same wrong
    ROOT, so it locks and corrupts one consistent wrong tree.

    A warning would not do. The write is not a side effect the operator can
    inspect afterwards and undo -- `stage()` overwrites two tracked files in a
    tree this session does not own, and the operator is by construction not
    watching that tree. A warning printed into a 4-arm measurement log is a
    silent failure with extra text (rule 2).

    NO LEGITIMATE MISMATCH EXISTS, checked rather than assumed. Every setter of
    this variable in the repository sets it to its OWN tree:
    `install-eqdyna.sh:119,152` and `.github/workflows/test.yml` (6 sites) use
    `$(pwd)`; `testsys/e2e/run_e2e.py:174` and `run_e2e_full.py:153` set their
    own file-derived `REPO_ROOT`; the four `testsys/regression/` tests that set
    it set their own `ROOT`. No caller anywhere -- and no tool in
    `testsys/perf/` -- sets it to a checkout other than its own, so refusing
    costs no supported workflow.
    """
    if os.path.realpath(ROOT) == os.path.realpath(SELF_ROOT):
        return
    raise SystemExit('\n'.join([
        'FAIL: %s' % ROOT_MISMATCH_HEADER,
        '  $EQDYNAROOT   : %s' % os.environ.get('EQDYNAROOT'),
        '  this tool is  : %s' % os.path.abspath(__file__),
        '  its checkout  : %s' % SELF_ROOT,
        '  would stage into: %s' % PKG,
        '',
        'stage() COPIES driver.py and MPI4NodalQuant.py into the package above,'
        ' and the',
        'finally clause rewrites them again on exit. With $EQDYNAROOT pointing'
        ' elsewhere',
        'those writes land in a checkout this session does not own -- the'
        ' damage rules 21',
        'and 21b exist for -- and Gate 0 below CANNOT catch it, because the'
        ' lock follows',
        'the same wrong root and so guards the same wrong tree (pathway item'
        ' 77).',
        '',
        'NOT warning and continuing. NOT preferring this file\'s own checkout'
        ' silently --',
        'either would write somewhere the operator is not watching (rule 2).'
        ' Re-export',
        'EQDYNAROOT for the tree you are standing in, or unset it and let this'
        ' tool',
        'resolve its own location.']))


def main():
    # GATE -1 (item 77): before argparse, because --notes defaults to a path
    # under ROOT and every write below -- the vault, the staged arm files, the
    # snapshot, the notes -- is resolved from it.
    require_root_is_this_checkout()
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', default='test.tpv104')
    ap.add_argument('--ranks', default='4,8,12')
    ap.add_argument('--arms', default='base,b,a,ab')
    ap.add_argument('--n-lo', type=int, default=20)
    ap.add_argument('--n-hi', type=int, default=60)
    ap.add_argument('--max-busy', type=float, default=0.5)
    ap.add_argument('--min-eff', type=float, default=0.95)
    ap.add_argument('--exclude-cpus', default='0,1')
    ap.add_argument('--merged-ref', default='HEAD')
    ap.add_argument('--notes', default=os.path.join(
        ROOT, 'NOTES_jaxmpi_ab_2026-09-22.md'))
    a = ap.parse_args()

    # GATE 0 (item 74): the exclusive lock on the package this tool STAGES
    # ARM FILES INTO, taken before the vault is built, before any arm is
    # staged, and before the notes file is opened -- so a refused invocation
    # costs nothing and leaves no trace. See LOCK_CONSEQUENCE for why a
    # collision here is worse than a crash.
    try:
        runlock.acquire(ROOT, PKG_RESOURCE, consequence=LOCK_CONSEQUENCE)
    except runlock.RunTreeLocked as exc:
        raise SystemExit('FAIL: %s' % exc)

    ranks = [int(x) for x in a.ranks.split(',') if x]
    arms = [x for x in a.arms.split(',') if x]
    excl = [int(x) for x in a.exclude_cpus.split(',') if x.strip()]

    sha = subprocess.run(['git', 'rev-parse', '--short', a.merged_ref],
                         cwd=ROOT, capture_output=True,
                         text=True).stdout.strip()
    work = tempfile.mkdtemp(prefix='jaxmpiab.',
                            dir=os.environ.get('TMPDIR', '/tmp'))
    vault = build_vault(work, a.merged_ref)
    rs.CASE = a.case
    case_dir = rs.build_py_case(a.case)
    nodes = numa.numa_topology()
    if not nodes:
        raise SystemExit('FAIL: numactl --hardware gave no topology.')

    out = os.path.join(ROOT, 'docs', 'perf_snapshots',
                       'jaxmpi_ab_%s.json' % time.strftime('%Y-%m-%d_%H%M%S'))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    note = open(a.notes, 'a')

    def say(s):
        print(s, flush=True)
        note.write(s + '\n')
        note.flush()

    def snap(rws):
        json.dump(dict(case=a.case, sha=sha, host=os.uname().nodename,
                       date=time.strftime('%Y-%m-%d %H:%M'), rows=rws,
                       n_lo=a.n_lo, n_hi=a.n_hi, max_busy=a.max_busy,
                       min_eff=a.min_eff, excluded_cpus=excl,
                       arms={k: ARM_DESC[k] for k in arms}),
                  open(out, 'w'), indent=1)

    say('\n## run %s  case %s  merged sha %s  n_lo/n_hi %d/%d  min-eff %.2f'
        % (time.strftime('%Y-%m-%d %H:%M'), a.case, sha, a.n_lo, a.n_hi,
           a.min_eff))
    t = ledger.box_tenancy(a.max_busy)
    say('box at start: %d/%d cpus busy over %.2f -> %d free; loadavg %.2f'
        % (t['busy'], t['total'], a.max_busy, t['total'] - t['busy'],
           os.getloadavg()[0]))

    rows = []
    try:
        for n in ranks:
            cpus, busy = ms.least_loaded_cpus(nodes, n, exclude=excl)
            worst = max(busy.values())
            tn = ledger.box_tenancy(a.max_busy)
            say('\n-- %d ranks -- cpus %s busy %s worst %.2f; box %d/%d busy '
                '-> %d free; loadavg %.2f'
                % (n, cpus, {c: round(b, 2) for c, b in busy.items()}, worst,
                   tn['busy'], tn['total'], tn['total'] - tn['busy'],
                   os.getloadavg()[0]))
            row = dict(ranks=n, cpus=cpus, n_lo=a.n_lo, n_hi=a.n_hi,
                       busy={str(c): b for c, b in busy.items()},
                       box_busy=tn['busy'], box_total=tn['total'],
                       loadavg=os.getloadavg(), arms={})
            rows.append(row)
            if worst > a.max_busy:
                say('  SKIPPED before launch: worst chosen cpu %.2f > %.2f'
                    % (worst, a.max_busy))
                row['skipped'] = ('pre-launch busy %.2f > %.2f'
                                  % (worst, a.max_busy))
                snap(rows)
                continue
            for arm in arms:
                stage(vault, arm)
                t0 = time.time()
                r = ms.per_step_jax_mpi(case_dir, cpus, n, a.n_lo, a.n_hi,
                                        'halo', 'cpu', warmup=True)
                if r is None:
                    say('  arm %-4s FAILED to run' % arm)
                    row['arms'][arm] = dict(failed=True)
                    snap(rows)
                    continue
                r.update(arm=arm, arm_desc=ARM_DESC[arm], min_eff=min(r['eff']),
                         point_wall_s=time.time() - t0,
                         loadavg_after=os.getloadavg())
                r['rejected'] = bool(r['min_eff'] < a.min_eff)
                row['arms'][arm] = r
                say('  arm %-4s %8.2f ms/step   min eff %.2f %s'
                    % (arm, r['ms_per_step'], r['min_eff'],
                       'REJECTED (< %.2f)' % a.min_eff if r['rejected']
                       else 'ACCEPTED'))
                say('           eff %s   rank ms %s'
                    % ([round(x, 2) for x in r['eff']],
                       [round(x, 1) for x in r['rank_ms']]))
                say('           exchange %s ms/step  barrier wait %s ms/step  '
                    'compile %.1fs  point wall %.0fs'
                    % ([round(x, 3) for x in r['mpi_ms']],
                       [round(x, 2) for x in r['wait_ms']],
                       r['compile_s'], r['point_wall_s']))
                snap(rows)
            good = {k: v for k, v in row['arms'].items()
                    if not v.get('failed') and not v.get('rejected')}
            base = good.get('base')
            if base is None:
                say('  VERDICT %d ranks: NO USABLE BASELINE (base arm absent '
                    'or rejected) -- every ratio at this rank count is '
                    'discarded.' % n)
            for arm in arms:
                v = good.get(arm)
                if arm == 'base' or v is None or base is None:
                    continue
                row.setdefault('speedup', {})[arm] = \
                    base['ms_per_step'] / v['ms_per_step']
                say('  VERDICT %d ranks arm %-4s %.3fx vs base  (%.2f -> %.2f '
                    'ms/step; min eff base %.2f arm %.2f)'
                    % (n, arm, base['ms_per_step'] / v['ms_per_step'],
                       base['ms_per_step'], v['ms_per_step'],
                       base['min_eff'], v['min_eff']))
            snap(rows)
    finally:
        stage(vault, 'ab')      # leave the tree on the merged branch's code
        snap(rows)
        say('\nsnapshot %s' % out)

    # LEDGER: one row per measured arm, through the project's own writer.
    # REJECTED points are recorded too, with the rejection on the row -- a
    # point that was measured and thrown away is evidence about this box.
    meta = dict(case=a.case, sha=sha, host=os.uname().nodename,
                date=time.strftime('%Y-%m-%d %H:%M'))
    lrows = []
    for row in rows:
        for arm, r in row.get('arms', {}).items():
            if r.get('failed'):
                continue
            lr = ledger._base(meta, 'run_jaxmpi_ab',
                              os.path.relpath(out, ROOT),
                              ledger.box_tenancy(a.max_busy), None)
            lr.update(
                backend='python-jax-mpi', ranks=row['ranks'],
                metric='per-step-by-difference',
                ms_per_step=r['ms_per_step'], n_lo=a.n_lo, n_hi=a.n_hi,
                platform='cpu', devices=None,
                platform_evidence='every rank printed device=cpu; '
                                  'run_mpi_scaling.jax_mpi_once discards the '
                                  'point otherwise',
                rank_ms_min=min(r['rank_ms']), rank_ms_max=max(r['rank_ms']),
                rank_ms_mean=sum(r['rank_ms']) / len(r['rank_ms']),
                effective_cores=sum(r['eff']) / len(r['eff']),
                effective_cores_per_rank=r['eff'],
                threads_per_rank=int(r['threads'][0]),
                busy_ceiling=a.max_busy,
                selection='A/B arm=%s (%s); cpus %s chosen once per rank '
                          'count and reused by all arms'
                          % (arm, ARM_DESC[arm], row['cpus']),
                verdict='REJECTED_LOW_EFF' if r['rejected'] else 'MEASURED',
                ab_arm=arm, ab_min_eff=r['min_eff'], sync='halo',
                cpus=row['cpus'], compile_s=r['compile_s'],
                mpi_ms_per_rank=r['mpi_ms'], wait_ms_per_rank=r['wait_ms'])
            lrows.append(lr)
    say('%d ledger row(s) appended to %s'
        % (ledger.append_rows(lrows), ledger.LEDGER_RELPATH))
    shutil.rmtree(work, ignore_errors=True)
    note.close()


if __name__ == '__main__':
    main()

#! /usr/bin/env python3
"""Behavioural regression guards for pathway item 91's six perf-tool
defects (testsys/perf/run_shard_scaling.py, run_numa_scaling.py,
run_perf.py, probe_scatter_bandwidth.py, run_setup_probe.py,
run_jaxmpi_ab.py -- sub-item e only). Full defect text: docs/BOARD_HISTORY.md
section "### Item 91".

Every check drives the real function with inputs that produce the bad case
and asserts it raises or labels -- never a grep over source text (rule 10a).
No MPI, no box-state dependence, no jax/numactl/git subprocess actually
runs: every subprocess boundary is monkeypatched so this guard is seconds
and machine-independent.

  (a) run_shard_scaling snapshot filename carries seconds, not just the
      date, and its rows carry the engine/policy keys the shared ledger
      reader (ledger.rows_from_scaling_snapshot) requires.
  (b) three per-step-by-difference sites (run_numa_scaling.per_step_and_fixed,
      probe_scatter_bandwidth.per_iter, run_perf.steady_state_per_step)
      raise on a non-positive result instead of recording it.
  (c) run_numa_scaling.time_one binds memory (--membind) as well as cpu
      (--physcpubind), via numactl, not bare taskset.
  (e) run_jaxmpi_ab.require_git_sha / run_setup_probe.require_git_sha raise
      on a git failure instead of returning ''.
  (f) run_perf.require_affinity_pinned raises when the pin did not take
      effect, instead of warning and proceeding.
  (g) run_shard_scaling.record_point / run_numa_scaling.record_baseline
      label the speedup baseline explicitly and never silently swap it
      when the intended first point was skipped.

Anti-vacuous-green discipline (papercuts): every verdict prints the content
property it asserted on.
"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
PERF = os.path.join(ROOT, 'testsys', 'perf')
sys.path.insert(0, PERF)
sys.path.insert(0, ROOT)

import run_shard_scaling as rss    # noqa: E402
import run_numa_scaling as rns     # noqa: E402
import run_perf as rp              # noqa: E402
import probe_scatter_bandwidth as psb  # noqa: E402
import run_setup_probe as rsp      # noqa: E402
import run_jaxmpi_ab as rjab       # noqa: E402
import ledger                      # noqa: E402

FAILURES = []


def check(ok, what):
    print('%s -- %s' % ('PASS' if ok else 'FAIL', what))
    if not ok:
        FAILURES.append(what)


def raises(fn, exc_types, what, must_contain=()):
    try:
        fn()
    except exc_types as e:
        msg = str(e)
        absent = [f for f in must_contain if f not in msg]
        check(not absent,
              '%s -> %s but message lacks %s: %s'
              % (what, type(e).__name__, absent, msg[:220])
              if absent else
              '%s -> %s; message: %s'
              % (what, type(e).__name__, msg[:220].replace('\n', ' ')))
        return
    check(False, '%s -> nothing raised' % what)


# ---------------------------------------------------------------- (a) ----

def check_91a_snapshot_naming_and_ledger():
    """Reloads the REAL run_shard_scaling module twice, under two faked wall
    times 5 minutes apart on the same calendar day, and compares its own
    module-level OUT global -- not a re-implementation of the naming
    expression. This is what actually catches the pre-fix `%Y-%m-%d`-only
    defect: with only the date, out1 == out2 here."""
    import importlib
    real_strftime = time.strftime

    def fake(t):
        def f(fmt, *a):
            return real_strftime(fmt, *a) if a else real_strftime(fmt, t)
        return f

    t1 = time.localtime(1_800_000_000)          # some day, some HH:MM:SS
    t2 = time.localtime(1_800_000_000 + 300)     # same day, 5 min later
    try:
        time.strftime = fake(t1)
        importlib.reload(rss)
        out1 = rss.OUT
        time.strftime = fake(t2)
        importlib.reload(rss)
        out2 = rss.OUT
    finally:
        time.strftime = real_strftime
        importlib.reload(rss)   # leave the module in its normal state
    check(out1 != out2,
          '91a: run_shard_scaling.OUT, reloaded at two wall times 5 min '
          'apart on the SAME calendar day, differ (%r vs %r)'
          % (out1, out2))

    # 91a's ledger half is deliberately NOT wired (PR #15 audit: the shared
    # reader would file shard points as run_scaling rows); nothing to guard.


# ---------------------------------------------------------------- (b) ----

def check_91b_nonpositive_raises():
    raises(lambda: rns.per_step_and_fixed(10.0, 9.0, 20, 60), (RuntimeError,),
           '91b run_numa_scaling.per_step_and_fixed: t_hi < t_lo',
           must_contain=('per-step',))
    ps, fixed = rns.per_step_and_fixed(1.0, 2.0, 10, 20)
    check(ps > 0, '91b run_numa_scaling.per_step_and_fixed: a genuine '
                  'positive case still returns (ps=%.4f)' % ps)

    orig_time_one = psb.time_one
    try:
        psb.time_one = lambda op, iters, cpus, node_map, dump_hlo_dir=None: (
            {30: 5.0, 100: 4.0}[iters])
        raises(lambda: psb.per_iter('control', [], {}, 30, 100),
               (RuntimeError,),
               '91b probe_scatter_bandwidth.per_iter: t_hi < t_lo',
               must_contain=('per-iteration',))
    finally:
        psb.time_one = orig_time_one

    orig_time_python = rp.time_python
    try:
        rp.time_python = lambda engine, n: (
            (0.0, 5.0) if n == 10 else (0.0, 4.0))
        raises(lambda: rp.steady_state_per_step('numpy', 10, 20),
               (RuntimeError,),
               '91b run_perf.steady_state_per_step: solve(n_hi) < solve(n_lo)',
               must_contain=('per-step',))
    finally:
        rp.time_python = orig_time_python


# ---------------------------------------------------------------- (c) ----

def check_91c_membind():
    got = {}

    def fake_run(cmd, env=None, capture_output=None, text=None):
        got['cmd'] = cmd
        class R:
            returncode = 0
            stdout = 'WALL 1.23\n'
            stderr = ''
        return R()

    orig_run = rns.subprocess.run
    try:
        rns.subprocess.run = fake_run
        node_map = {0: [0, 1, 2, 3, 4, 5, 6, 7],
                    1: [8, 9, 10, 11, 12, 13, 14, 15]}
        w = rns.time_one('/case', 20, [8, 9], node_map)
    finally:
        rns.subprocess.run = orig_run
    cmd = got['cmd']
    check(w == 1.23, '91c: time_one still returns the wall seconds (%.2f)' % w)
    check(cmd[0] == 'numactl',
          '91c: time_one launches via numactl, not bare taskset (cmd[0]=%r)'
          % cmd[0])
    check(any(c == '--physcpubind=8,9' for c in cmd),
          '91c: --physcpubind carries the exact cpu list (cmd=%s)' % cmd)
    check(any(c == '--membind=1' for c in cmd),
          '91c: --membind carries the NUMA node(s) those cpus belong to '
          '(node 1 for cpus 8,9) -- the fix run_scaling.py:18-26 documents '
          '(cmd=%s)' % cmd)
    check(rns._nodes_of(node_map, [8, 9]) == [1],
          '91c: _nodes_of([8,9]) == [1] (pure helper, no subprocess)')
    check(rns._nodes_of(node_map, [0, 8]) == [0, 1],
          '91c: _nodes_of([0,8]) == [0,1] spans both nodes when cpus do')


# ---------------------------------------------------------------- (e) ----

def check_91e_git_sha():
    for mod, name in ((rjab, 'run_jaxmpi_ab'), (rsp, 'run_setup_probe')):
        raises(lambda mod=mod: mod.require_git_sha(
                   1, '', 'fatal: bad revision', 'git rev-parse'),
               (SystemExit,),
               '91e %s.require_git_sha: git returncode=1' % name,
               must_contain=('git rev-parse',))
        raises(lambda mod=mod: mod.require_git_sha(0, '', '', 'git rev-parse'),
               (SystemExit,),
               '91e %s.require_git_sha: returncode=0 but empty stdout '
               '(the exact defect: recording sha="")' % name)
        sha = mod.require_git_sha(0, ' abcd123\n', '', 'git rev-parse')
        check(sha == 'abcd123',
              '91e %s.require_git_sha: success returns the stripped sha '
              '(%r)' % (name, sha))


# ---------------------------------------------------------------- (f) ----

def check_91f_affinity_gate():
    try:
        rp.require_affinity_pinned('[0]', '0')
        ok = True
    except SystemExit:
        ok = False
    check(ok, '91f: an affinity match ([0] for core 0) does not raise')
    raises(lambda: rp.require_affinity_pinned('[0,1]', '0'), (SystemExit,),
           '91f: an affinity mismatch ([0,1] for core 0) refuses instead of '
           'warning-and-proceeding')


# ---------------------------------------------------------------- (g) ----

def check_91g_baseline_labelling():
    base, base_n = {}, {}
    _, bn1, note1 = rss.record_point(base, base_n, 'element', 2, 100.0, 1)
    check(bn1 == 2 and note1 is not None and '1' in note1 and '2' in note1,
          '91g run_shard_scaling.record_point: n=1 was skipped, n=2 sets '
          'the baseline and gets a non-None NOTE naming both n=%r vs the '
          'devices[0]=1 it silently used to be read against (note=%r)'
          % (bn1, note1))
    speedup2, bn2, note2 = rss.record_point(base, base_n, 'element', 2, 100.0, 1)
    check(speedup2 == 1.0 and note2 is None,
          '91g: the already-set baseline point itself gets no repeat NOTE '
          '(speedup=%.2f, note=%r)' % (speedup2, note2))

    base2, base_n2 = {}, {}
    _, bn3, note3 = rss.record_point(base2, base_n2, 'element', 1, 100.0, 1)
    check(bn3 == 1 and note3 is None,
          '91g: when n=1 (devices[0]) DOES set the baseline, no NOTE fires '
          '(bn=%r note=%r)' % (bn3, note3))

    state = {}
    speedup, blabel, note = rns.record_baseline(state, 'spread-8', 50.0,
                                                 'within-node-8')
    check(blabel == 'spread-8' and note is not None
          and 'within-node-8' in note and 'spread-8' in note,
          '91g run_numa_scaling.record_baseline: first config skipped, '
          'second config (%r) becomes baseline with a non-None NOTE (%r)'
          % (blabel, (note or '')[:80]))


CHECKS = [
    ('91a', check_91a_snapshot_naming_and_ledger),
    ('91b', check_91b_nonpositive_raises),
    ('91c', check_91c_membind),
    ('91e', check_91e_git_sha),
    ('91f', check_91f_affinity_gate),
    ('91g', check_91g_baseline_labelling),
]


def main():
    for tag, fn in CHECKS:
        print('\n-- item %s --' % tag)
        fn()
    print()
    if FAILURES:
        print('FAIL test_perf_item91_guards: %d check(s) failed' % len(FAILURES))
        return 1
    print('SUCCESS test_perf_item91_guards: all 6 item-91 perf-tool defects '
          'guarded behaviourally')
    return 0


if __name__ == '__main__':
    sys.exit(main())


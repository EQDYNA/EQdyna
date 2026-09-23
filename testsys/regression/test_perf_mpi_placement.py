#! /usr/bin/env python3
"""Guard for `--placement` on run_mpi_scaling.py and the `placement`
discriminator on run_mpi_scaling ledger rows (2026-09-23, mission: perf-
tooling guard).

Incident guarded: `least_loaded_cpus` sorts by (busy, node, cpu); on a quiet
box (most cpus read busy 0.00) the node tiebreak dominates and it PACKS ranks
onto the fewest NUMA nodes rather than spreading them. Measured, same
session, test.tpv104 16 ranks (docs/perf_ledger.jsonl lines 361-364): SPREAD
(2 ranks/node) jax-MPI 52.5 vs Fortran 63.3 ms/step = 0.83x; PACKED (6/7/3 on
3 nodes) 96.0 vs 56.0 = 1.71x. A whole day of 8/16-rank tables was
placement-biased with nothing on the row saying so.

Checks:
  1. topology behaviour (not source text, rule 10a) -- an 8-node x 8-cpu
     synthetic busy map, all cpus idle: packed puts k=16 ranks on 2 nodes;
     spread puts 2 on each of 8.
  2. one node fully busy -- spread excludes it and uses the remaining 7,
     ranks-per-node differing by at most 1.
  3. spread refuses (SystemExit) when too few cpus qualify under the ceiling.
  4. argparse: neither --placement nor --cpus given exits non-zero, and
     names the confound rather than defaulting silently.
  5. ledger: a run_mpi_scaling row missing `placement` is refused on append;
     the same row WITH placement (and ranks_per_node) appends cleanly.
  6. the converter (rows_from_mpi_scaling_snapshot) propagates placement and
     ranks_per_node from a synthetic snapshot row into the ledger row.

Anti-vacuous-green discipline (papercuts): every verdict prints the content
property it asserted on -- node counts, chosen cpu sets, exit codes.
"""
import copy
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
PERF = os.path.join(ROOT, 'testsys', 'perf')
sys.path.insert(0, PERF)
import ledger              # noqa: E402
import run_mpi_scaling as rms  # noqa: E402

FAILURES = []


def check(ok, what):
    print('%s -- %s' % ('PASS' if ok else 'FAIL', what))
    if not ok:
        FAILURES.append(what)


def refuses(fn, what, must_contain=()):
    try:
        fn()
    except SystemExit as e:
        msg = str(e)
        absent = [f for f in must_contain if f not in msg]
        check(not absent,
              '%s -> SystemExit%s; message: %s'
              % (what, '' if not absent else ' but message lacks %s' % absent,
                 msg[:200].replace('\n', ' ')))
        return
    check(False, '%s -> nothing raised' % what)


# Synthetic 8-node x 8-cpu topology (this box happens to be exactly this
# shape too, but the test does not rely on that -- it builds its own map).
NODES = {n: list(range(n * 8, n * 8 + 8)) for n in range(8)}


def counts_by_node(cpus):
    cpu2node = {c: n for n, cs in NODES.items() for c in cs}
    out = {}
    for c in cpus:
        out[cpu2node[c]] = out.get(cpu2node[c], 0) + 1
    return out


def main():
    # 1. quiet box: packed packs, spread spreads.
    quiet = {c: 0.0 for cs in NODES.values() for c in cs}
    packed = rms._pick_packed(NODES, 16, quiet)
    spread = rms._pick_spread(NODES, 16, quiet, max_busy=0.5)
    pb, sb = counts_by_node(packed), counts_by_node(spread)
    check(len(pb) == 2 and sorted(pb.values()) == [8, 8],
          'packed(16, quiet box): per-node counts %r (2 nodes, 8 each)' % pb)
    check(len(sb) == 8 and sorted(sb.values()) == [2] * 8,
          'spread(16, quiet box): per-node counts %r (8 nodes, 2 each)' % sb)

    # 2. one node fully busy: spread excludes it, uses the remaining 7,
    # counts differ by at most 1.
    busy_map = dict(quiet)
    for c in NODES[3]:
        busy_map[c] = 1.0
    spread2 = rms._pick_spread(NODES, 16, busy_map, max_busy=0.5)
    b2 = counts_by_node(spread2)
    check(3 not in b2 and len(b2) == 7
          and max(b2.values()) - min(b2.values()) <= 1,
          'spread(16, node 3 fully busy): per-node counts %r -- node 3 '
          'excluded, 7 nodes used, counts differ by <=1' % b2)

    # 3. spread refuses when too few cpus qualify.
    scarce = dict(quiet)
    for n in (1, 2, 3, 4, 5, 6, 7):
        for c in NODES[n]:
            scarce[c] = 1.0
    refuses(lambda: rms._pick_spread(NODES, 16, scarce, max_busy=0.5),
            'spread(16, 7 of 8 nodes busy, only 8 cpus qualify)',
            must_contain=('8', '16'))

    # 4. argparse: neither --placement nor --cpus -> non-zero exit, names
    # the confound. Invoked as a subprocess so the early SystemExit is
    # observed before any heavy work (build_py_case, mpirun) could run --
    # this is the fast fail-closed path, not a live measurement.
    r = subprocess.run(
        [sys.executable, os.path.join(PERF, 'run_mpi_scaling.py'),
         '--ranks', '1', '--skip-fortran'],
        capture_output=True, text=True, cwd=ROOT,
        env=dict(os.environ, EQDYNAROOT=ROOT))
    check(r.returncode != 0
          and 'placement' in r.stderr and 'PACKS' in r.stderr,
          'CLI with neither --placement nor --cpus: exit %d, stderr names '
          'the confound (%s)' % (r.returncode, r.stderr.strip()[-200:]))

    # 5. ledger: placement required on append for run_mpi_scaling rows.
    tmp = tempfile.mkdtemp(prefix='perf_mpi_placement_test.')
    led = os.path.join(tmp, 'ledger.jsonl')
    base_row = dict(
        ts_utc='2026-09-23T12:00:00Z', snapshot_date_local='2026-09-23 12:00',
        sha='0123abc', host='cotopaxi', tool='run_mpi_scaling',
        case='test.tpv104', backend='python-jax-mpi', ranks=16,
        ms_per_step=52.5, rank_ms_min=50.1, rank_ms_max=55.2,
        rank_ms_mean=52.4, effective_cores=3.9, threads_per_rank=1,
        busy_ceiling=0.6, tenancy_busy=20, tenancy_total=64,
        metric='per-step-by-difference', n_lo=40, n_hi=160,
        snapshot='docs/perf_snapshots/mpi_scaling_2026-09-23_112435_tpv104.json',
        platform='cpu', devices=None,
        platform_evidence='requested cpu: JAX_PLATFORMS=cpu pins the jax '
                          'backend to host devices',
        parallelism='mpi')
    missing = copy.deepcopy(base_row)
    refuses_val = None
    try:
        ledger.append(missing, path=led)
        check(False, 'append of a run_mpi_scaling row WITHOUT placement')
    except ValueError as e:
        refuses_val = str(e)
        check('placement' in refuses_val,
              'append of a run_mpi_scaling row WITHOUT placement -> '
              'ValueError: %s' % refuses_val[:200])
    nlines = sum(1 for _ in open(led)) if os.path.exists(led) else 0
    check(nlines == 0,
          'no refused row reached the ledger (%d lines)' % nlines)

    with_placement = dict(base_row, placement='spread',
                          ranks_per_node=[2, 2, 2, 2, 2, 2, 2, 2])
    ledger.append(with_placement, path=led)
    got = [json.loads(l) for l in open(led)]
    check(len(got) == 1 and got[0]['placement'] == 'spread'
          and got[0]['ranks_per_node'] == [2, 2, 2, 2, 2, 2, 2, 2],
          'the same row WITH placement appends: 1 line, placement=%r '
          'ranks_per_node=%r' % (got[0]['placement'], got[0]['ranks_per_node']))

    # a row from a DIFFERENT tool never needs placement.
    other_tool = dict(base_row, tool='run_scaling')
    if 'placement' in other_tool:
        del other_tool['placement']
    ledger.append(other_tool, path=led)
    check(sum(1 for _ in open(led)) == 2,
          'a run_scaling row with no placement key appends fine (2 lines) '
          '-- the field is scoped to run_mpi_scaling only')

    # a legacy run_mpi_scaling row (predates this field) still validates on
    # READ, never on append.
    legacy = copy.deepcopy(base_row)
    try:
        ledger.validate(legacy)
        check(True, 'legacy run_mpi_scaling row (no placement key) validates '
                    'when READ')
    except ValueError as e:
        check(False, 'legacy run_mpi_scaling row validates when READ (%s)' % e)

    # 6. converter propagation from a synthetic snapshot.
    mpi_meta = dict(case='test.tpv104', sha='abc1234', host='cotopaxi',
                    date='2026-09-23 16:37', n_lo=40, n_hi=160, max_busy=0.6,
                    rows=[dict(ranks=16, cpus=list(range(16)), n_lo=40, n_hi=160,
                               placement='spread',
                               ranks_per_node=[2, 2, 2, 2, 2, 2, 2, 2],
                               jax_halo=dict(
                                   ms_per_step=52.5,
                                   rank_ms=[50.1, 55.2, 52.0, 52.5] * 4,
                                   eff=[3.9] * 16, threads=[1] * 16),
                               fortran=dict(ms_per_step=63.3))])
    conv = ledger.rows_from_mpi_scaling_snapshot(
        mpi_meta, 'docs/perf_snapshots/mpi_scaling_x.json',
        dict(busy=20, total=64))
    check(len(conv) == 2
          and all(r['placement'] == 'spread' for r in conv)
          and all(r['ranks_per_node'] == [2, 2, 2, 2, 2, 2, 2, 2]
                  for r in conv),
          'converter propagates placement=%s ranks_per_node=%s onto both '
          'the jax-mpi and the fortran row'
          % ({r['backend']: r['placement'] for r in conv},
             {r['backend']: r['ranks_per_node'] for r in conv}))
    for r in conv:
        ledger.append(r, path=led)
    check(sum(1 for _ in open(led)) == 4,
          'both converted rows append cleanly (4 lines total)')

    print()
    if FAILURES:
        print('FAIL test_perf_mpi_placement: %d check(s) failed' % len(FAILURES))
        return 1
    print('SUCCESS test_perf_mpi_placement: packed/spread selection is '
          'behaviour-tested and every run_mpi_scaling ledger row must say '
          'which one produced it')
    return 0


if __name__ == '__main__':
    sys.exit(main())

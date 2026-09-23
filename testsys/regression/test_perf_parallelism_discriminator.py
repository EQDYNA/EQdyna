#! /usr/bin/env python3
"""Guard for the `parallelism` discriminator on perf-ledger rows
(testsys/perf/ledger.py, docs/perf_ledger.jsonl).

Incident guarded (2026-09-23, found by the owner reading the board table):
`ranks` carried two incommensurable quantities in ONE column with nothing on
the row saying which. On a run_scaling `python-jax` row it is the cpus in ONE
process's affinity mask (run_scaling.time_one_py launches a single
`sys.executable -c` under numactl; XLA and OpenBLAS auto-size their pools from
the pin). On a run_scaling `fortran` row, and on every run_mpi_scaling and
run_jaxmpi_ab row, it is `mpirun -np n` -- n separate processes. Dividing one
by the other published "Fortran 16.24x vs jax 2.59x at 16"; the real-ranks
measurement put jax at 1.58x of Fortran, not 4.19x.

The fix is a REQUIRED, non-defaulted field, and this guard pins it:
  1. the value set is exactly ('mpi', 'threads', 'serial');
  2. appending a row WITHOUT the field is refused, with a message that names
     the emitter, and nothing reaches the ledger;
  3. an out-of-set value, and a 'serial' row claiming ranks > 1, are refused;
  4. a legacy row (no key) still validates when READ -- the ledger is
     append-only, history is never rewritten to gain the field, and absence
     means UNKNOWN, never a default;
  5. per emitter, the value is the right one for how that emitter launches:
     run_scaling fortran=mpi / python-*=threads, run_mpi_scaling both=mpi,
     run_e2e mpirun cells=mpi and single-process cells=serial;
  6. an engine/backend the converters have not been taught RAISES instead of
     being labelled by guess;
  7. every emitter in the tree is covered -- a NEW caller of append_rows that
     does not set the field fails this check before it fails at runtime;
  8. the committed ledger: old rows keyless (UNKNOWN), new rows legal.

Anti-vacuous-green discipline (papercuts): every verdict prints the content
property it asserted on -- the parsed value, the refusal message, the counts.
"""
import copy
import json
import os
import re
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
PERF = os.path.join(ROOT, 'testsys', 'perf')
sys.path.insert(0, PERF)
import ledger  # noqa: E402

FAILURES = []


def check(ok, what):
    print('%s -- %s' % ('PASS' if ok else 'FAIL', what))
    if not ok:
        FAILURES.append(what)


def refuses(fn, what, must_contain=()):
    """Call fn; PASS only if it raises ValueError whose message contains every
    fragment in must_contain. Prints the message, so a refusal that fires for
    an unrelated reason is visible rather than counted as a green."""
    try:
        fn()
    except ValueError as e:
        msg = str(e)
        absent = [f for f in must_contain if f not in msg]
        check(not absent,
              '%s -> ValueError%s; message: %s'
              % (what, '' if not absent else ' but message lacks %s' % absent,
                 msg[:200].replace('\n', ' ')))
        return
    check(False, '%s -> nothing raised (the field is being defaulted in)'
                 % what)


def good_row(**kw):
    """A row shaped like what an emitter produces, with the discriminator."""
    r = dict(ts_utc='2026-09-23T12:00:00Z',
             snapshot_date_local='2026-09-23 12:00',
             sha='0123abc', host='cotopaxi', tool='run_scaling',
             case='test.tpv104', backend='fortran', ranks=16,
             ms_per_step=64.74, rank_ms_min=None, rank_ms_max=None,
             rank_ms_mean=None, effective_cores=None, threads_per_rank=None,
             busy_ceiling=0.5, tenancy_busy=20, tenancy_total=64,
             metric='per-step-by-difference', n_lo=20, n_hi=60,
             snapshot='docs/perf_snapshots/scaling_2026-09-23_x.json',
             platform='cpu', devices=None,
             platform_evidence='synthetic guard row',
             parallelism='mpi')
    r.update(kw)
    return r


SCALING_META = dict(
    case='test.tpv104', sha='abc1234', host='cotopaxi',
    date='2026-09-23 11:35', busy_ceiling=0.45,
    rows=[dict(engine='fortran', n=16, policy='compact', ms_per_step=64.74,
               n_lo=20, n_hi=60, cpus=list(range(16))),
          dict(engine='python-jax', n=16, policy='compact', ms_per_step=271.5,
               n_lo=20, n_hi=60, cpus=list(range(16))),
          dict(engine='python-numpy', n=8, policy='spread', ms_per_step=250.25,
               n_lo=20, n_hi=60, cpus=[0, 8, 16, 24, 32, 40, 48, 56])])

MPI_META = dict(
    case='test.tpv104', sha='abc1234', host='cotopaxi',
    date='2026-09-23 16:37', n_lo=40, n_hi=160, max_busy=0.6,
    rows=[dict(ranks=4, cpus=[8, 9, 16, 17], n_lo=40, n_hi=160,
               jax_halo=dict(ms_per_step=605.0,
                             rank_ms=[600.0, 605.0, 590.0, 601.0],
                             eff=[3.9, 3.8, 4.0, 3.7], threads=[1, 1, 1, 1]),
               fortran=dict(ms_per_step=101.5))])

E2E_META = dict(
    tool='run_e2e', sha='abc1234', host='cotopaxi', date='2026-09-23 07:00',
    label='default', device='cpu', jobs_budget=8, tenancy_ceiling=0.5,
    cells=[dict(case='test.tpv8', backend='fortran', ok=True, seconds=41.2,
                ranks=4, platform='cpu',
                platform_evidence='fortran: no GPU path'),
           dict(case='test.tpv8', backend='python-jax-mpi', ok=True,
                seconds=55.3, ranks=4, platform='cpu',
                platform_evidence='every rank printed device=cpu'),
           dict(case='test.tpv8', backend='python-numpy', ok=True,
                seconds=97.6, ranks=1, platform='cpu',
                platform_evidence='numpy: no GPU path'),
           dict(case='test.tpv8', backend='python-jax', ok=True, seconds=63.0,
                ranks=1, platform='cpu',
                platform_evidence='JAX_PLATFORMS=cpu pinned')])


def main():
    tmp = tempfile.mkdtemp(prefix='perf_parallelism_test.')
    led = os.path.join(tmp, 'ledger.jsonl')

    # 1. the value set, pinned. Widening it is a deliberate edit that must
    # also update the reader contract in the module docstring.
    check(ledger.PARALLELISM == ('mpi', 'threads', 'serial'),
          'value set is exactly %r' % (ledger.PARALLELISM,))

    # 2. fail closed on a missing field, naming the emitter.
    absent = good_row()
    del absent['parallelism']
    refuses(lambda: ledger.append(absent, path=led),
            'append of a row WITHOUT parallelism (tool=run_scaling)',
            must_contain=("'parallelism'", "tool='run_scaling'",
                          "case='test.tpv104'"))

    # 3. out-of-set value, and the serial/ranks contradiction.
    refuses(lambda: ledger.append(good_row(parallelism='openmp'), path=led),
            "append with parallelism='openmp'",
            must_contain=('openmp', 'mpi'))
    refuses(lambda: ledger.append(good_row(parallelism=None), path=led),
            'append with parallelism=None (a null is not an opt-out)',
            must_contain=('None',))
    refuses(lambda: ledger.append(
                good_row(parallelism='serial', ranks=16), path=led),
            "append with parallelism='serial' but ranks=16",
            must_contain=('serial', '16'))
    nlines = sum(1 for _ in open(led)) if os.path.exists(led) else 0
    check(nlines == 0,
          'no refused row reached the ledger (%d lines in the temp ledger)'
          % nlines)

    # ... and the well-formed row of the same shape DOES land, so the refusals
    # above are about the field and not about the fixture being unappendable.
    ledger.append(good_row(), path=led)
    got = [json.loads(l) for l in open(led)]
    check(len(got) == 1 and got[0]['parallelism'] == 'mpi',
          'the same row WITH parallelism appends: 1 line, parallelism=%r, '
          'ranks=%r' % (got[0]['parallelism'], got[0]['ranks']))

    # 4. legacy rows: absent key is UNKNOWN when READ, never defaulted.
    legacy = good_row()
    del legacy['parallelism']
    try:
        ledger.validate(legacy)
        check('parallelism' not in legacy,
              'a legacy row (no parallelism key) validates when READ and '
              'gains no key -- absence is UNKNOWN, not "mpi"')
    except ValueError as e:
        check(False, 'legacy row validates when READ (%s)' % e)

    # 5. per-emitter values, from each converter on recorded-shape input.
    srows = ledger.rows_from_scaling_snapshot(
        SCALING_META, 'docs/perf_snapshots/scaling_x.json',
        dict(busy=10, total=64))
    smap = {r['backend']: r['parallelism'] for r in srows}
    check(smap == {'fortran': 'mpi', 'python-jax': 'threads',
                   'python-numpy': 'threads'},
          'emitter run_scaling: %r -- fortran is mpirun -np n, the python '
          'engines are ONE process under a numactl mask (the defect)' % smap)
    jax16 = next(r for r in srows if r['backend'] == 'python-jax')
    f16 = next(r for r in srows if r['backend'] == 'fortran')
    check(jax16['ranks'] == f16['ranks'] == 16
          and jax16['parallelism'] != f16['parallelism'],
          'emitter run_scaling: the two rows that were divided into the bogus '
          '4.19x both say ranks=16 but now say parallelism=%r vs %r, so a '
          'reader can see they are not the same quantity'
          % (f16['parallelism'], jax16['parallelism']))

    mrows = ledger.rows_from_mpi_scaling_snapshot(
        MPI_META, 'docs/perf_snapshots/mpi_scaling_x.json',
        dict(busy=30, total=64))
    mmap = {r['backend']: r['parallelism'] for r in mrows}
    check(len(mrows) == 2 and mmap == {'python-jax-mpi': 'mpi',
                                       'fortran': 'mpi'},
          'emitter run_mpi_scaling: %d row(s), %r -- both backends are '
          'mpirun -np ranks' % (len(mrows), mmap))

    erows = ledger.rows_from_e2e_results(
        E2E_META, 'docs/perf_snapshots/e2e_cells_x.json',
        dict(busy=12, total=64))
    # `.get`, not `[...]`: an emitter that stops setting the field must fail
    # with the verdict below naming the emitter and the row, not with a bare
    # KeyError traceback from inside this comprehension. A guard that cannot
    # say what it found teaches the reader to re-run instead of read.
    emap = {r['backend']: (r.get('parallelism', '<MISSING>'), r['ranks'])
            for r in erows}
    check(emap == {'fortran': ('mpi', 4), 'python-jax-mpi': ('mpi', 4),
                   'python-numpy': ('serial', 1), 'python-jax': ('serial', 1)},
          'emitter run_e2e: %r -- the mpirun cells are mpi, the '
          'single-process cells are serial with ranks=1 counting the PROCESS '
          '(not one core)' % emap)

    # 6. an unknown engine/backend is refused, not labelled by guess.
    unknown_engine = copy.deepcopy(SCALING_META)
    unknown_engine['rows'] = [dict(SCALING_META['rows'][0],
                                   engine='python-torch')]
    refuses(lambda: ledger.rows_from_scaling_snapshot(
                unknown_engine, 'docs/perf_snapshots/scaling_x.json',
                dict(busy=10, total=64)),
            "run_scaling converter on an undeclared engine 'python-torch'",
            must_contain=('python-torch',))
    unknown_cell = copy.deepcopy(E2E_META)
    unknown_cell['cells'] = [dict(E2E_META['cells'][0],
                                  backend='python-torch')]
    refuses(lambda: ledger.rows_from_e2e_results(
                unknown_cell, 'docs/perf_snapshots/e2e_cells_x.json',
                dict(busy=12, total=64)),
            "run_e2e converter on an undeclared backend 'python-torch'",
            must_contain=('python-torch',))

    # 7. emitter coverage: every file that appends must set the field, either
    # through one of the guarded converters or with a literal of its own.
    # A new tool added without it fails HERE, in seconds, not mid-measurement.
    # `capture_e2e_cells[_or_warn]` counts as covered: it is ledger.py's own
    # wrapper around the guarded rows_from_e2e_results converter.
    CONVERTERS = ('rows_from_scaling_snapshot',
                  'rows_from_mpi_scaling_snapshot', 'rows_from_e2e_results',
                  'capture_e2e_cells')
    appenders, uncovered = [], []
    for dirpath, dirnames, filenames in os.walk(os.path.join(ROOT, 'testsys')):
        dirnames[:] = [d for d in dirnames if d != '__pycache__']
        for fn in sorted(filenames):
            if not fn.endswith('.py'):
                continue
            p = os.path.join(dirpath, fn)
            rel = os.path.relpath(p, ROOT)
            if rel.startswith(os.path.join('testsys', 'regression')):
                continue            # guards, not emitters
            if rel == os.path.join('testsys', 'perf', 'ledger.py'):
                continue            # the module under guard itself
            src = open(p).read()
            if not re.search(r'\bledger\.(append(_rows)?|capture_e2e_cells'
                             r'(_or_warn)?)\(', src):
                continue
            appenders.append(rel)
            via_converter = any(c in src for c in CONVERTERS)
            literal = re.search(r"parallelism=['\"](%s)['\"]"
                                % '|'.join(ledger.PARALLELISM), src)
            if not (via_converter or literal):
                uncovered.append(rel)
    check(bool(appenders) and not uncovered,
          'emitter coverage: %d appending tool(s) %s; each goes through a '
          'guarded converter or sets a literal parallelism= (uncovered: %s)'
          % (len(appenders), appenders, uncovered or 'none'))
    ab = os.path.join(ROOT, 'testsys', 'perf', 'run_jaxmpi_ab.py')
    ab_src = open(ab).read()
    check("parallelism='mpi'" in ab_src,
          'emitter run_jaxmpi_ab (builds its row inline, no converter): '
          "source sets parallelism='mpi' -- every arm is mpirun -np ranks")

    # 7b. mutation proof, in-suite: strip the field from a REAL converter row
    # and the append must refuse; the unstripped row appends.
    real = copy.deepcopy(jax16)
    del real['parallelism']
    refuses(lambda: ledger.append(real, path=led),
            'append of a real run_scaling python row with the field stripped',
            must_contain=("'parallelism'", "tool='run_scaling'"))
    before = sum(1 for _ in open(led))
    ledger.append(copy.deepcopy(jax16), path=led)
    after = [json.loads(l) for l in open(led)]
    check(len(after) == before + 1 and after[-1]['parallelism'] == 'threads',
          'the same converter row unmodified appends (%d -> %d lines, last '
          'parallelism=%r)' % (before, len(after), after[-1]['parallelism']))

    # 8. the committed ledger: old rows keyless, any present value legal.
    real_led = os.path.join(ROOT, 'docs', 'perf_ledger.jsonl')
    if not os.path.exists(real_led):
        check(False, 'committed ledger docs/perf_ledger.jsonl exists')
    else:
        rows = [json.loads(l) for l in open(real_led) if l.strip()]
        withf = [r for r in rows if 'parallelism' in r]
        without = [r for r in rows if 'parallelism' not in r]
        bad = [(i + 1, r.get('parallelism')) for i, r in enumerate(rows)
               if 'parallelism' in r
               and r['parallelism'] not in ledger.PARALLELISM]
        check(not bad,
              'committed ledger: %d row(s), %d carry parallelism, %d are '
              'legacy/UNKNOWN; illegal values: %s'
              % (len(rows), len(withf), len(without), bad or 'none'))
        unreadable = []
        for i, r in enumerate(rows, 1):
            try:
                ledger.validate(r)
            except ValueError as e:
                unreadable.append('line %d: %s' % (i, e))
        check(not unreadable,
              'committed ledger: all %d row(s) still validate when READ -- '
              'legacy rows are not required to gain the field (problems: %s)'
              % (len(rows), unreadable or 'none'))

    print()
    if FAILURES:
        print('FAIL test_perf_parallelism_discriminator: %d check(s) failed'
              % len(FAILURES))
        return 1
    print('SUCCESS test_perf_parallelism_discriminator: `ranks` now says what '
          'it counts, on every emitter, fail-closed')
    return 0


if __name__ == '__main__':
    sys.exit(main())

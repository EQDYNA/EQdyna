#! /usr/bin/env python3
"""
THE test. One sweep, two nested loops, one comparison:

    for case in CASES:                       # testNameList.py, 8 cases
        for backend in BACKENDS:             # fortran, python-numpy, python-jax
            run(case, backend) -> canonical frt -> compare vs the ONE
                                  committed reference, at THAT CASE's bound

`backend` is an axis of the sweep, not a tier. The former `accept` tier is
gone: it ran the SAME solver against the SAME references with a SECOND
comparison implementation, over 5 listed cases while this tier gated 8, on the
jax backend only because it invoked the standalone with no --backend. Its
"SUCCESS accept (5/5 cases)" therefore reported five of five LISTED -- 5 of 16
case x backend combinations -- and read like completeness. Both of its halves
live here now: the run half below, the comparison half in testsys/compare.py.

WHAT A GREEN RUN MEANS, EXACTLY
It means every cell this run PRINTED as "will RUN" was run and gated. The
coverage block is printed before the first case starts and restated in the
summary, so the result can never be read as broader than it is. Cells outside
the selection are named, and cells the table declares unsupported are named
with their reason.

DECLARED-UNSUPPORTED IS NOT A SKIP
  - default sweep (no --cases/--backends): the table's unsupported cells are
    printed with their reasons and are NOT counted as coverage.
  - explicit selection naming an unsupported cell: hard FAILURE. The caller
    expected that cell to run; "I could not check this" and "this is fine"
    must not share an exit code (PROJECT_RULES rule 2).

Named gates, in order:
  1. testsys/regression/test_create_newcase.py -- cheap check first (rule 9).
  2. Fresh build via ./install-eqdyna.sh (only when a fortran cell is selected;
     the previous bin/eqdyna is removed first so a stale binary cannot paper
     over a broken build). EQDYNA_E2E_BIN overrides with a pre-built binary
     when a concurrent long job must not have bin/ disturbed -- it skips the
     build-and-install path, so a green run under it says nothing about that
     path (rule 16); it is printed loudly.
  3. Every selected cell: run, then compare.

Rule 8: the previous test/ tree is rotated to test.prev/ rather than deleted,
and every cell's output stays on disk under test/ -- including the python
cells, which the old accept tier ran in a tempfile directory and deleted, so a
failure destroyed its own evidence.

No timeouts: this is a shared box and dynamic-rupture runs are slow when it is
busy.
"""
import argparse
import os
import shutil
import subprocess
import sys
import time

E2E_DIR = os.path.dirname(os.path.abspath(__file__))
TESTSYS = os.path.dirname(E2E_DIR)
REPO_ROOT = os.path.dirname(TESTSYS)
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from testsys import compare, matrix  # noqa: E402

# Line-buffered stdout. Redirected to a file or through `tee`, Python block-
# buffers its OWN prints while subprocess children write straight to the fd --
# so the coverage block, the per-cell verdicts and the child output interleave
# wrongly, and none of this script's own lines survive a kill. A GitHub runner
# SIGTERM'd one of these runs (exit 143) and the log showed a silence and an
# exit code. The block that states what a run covered is the last thing that
# should be lost when a run dies.
sys.stdout.reconfigure(line_buffering=True)

MACHINE = os.environ.get('EQDYNA_TEST_MACHINE', 'ubuntu')
MPIRUN = os.environ.get('EQDYNA_MPIRUN', 'mpirun')
BIN_OVERRIDE = os.environ.get('EQDYNA_E2E_BIN')


def _run(cmd, cwd, env):
    print('+ (%s) %s' % (os.path.basename(cwd), ' '.join(cmd)))
    return subprocess.call(cmd, cwd=cwd, env=env)


def base_env():
    env = dict(os.environ)
    env['EQDYNAROOT'] = REPO_ROOT
    env['PATH'] = os.pathsep.join([
        os.path.join(REPO_ROOT, 'bin'),
        os.path.join(REPO_ROOT, 'scripts'),
        env.get('PATH', ''),
    ])
    return env


# --------------------------------------------------------------------------
# cell runners
# --------------------------------------------------------------------------
def make_serial_case(case_name, case_dir, env):
    """create.newcase + force a SERIAL decomposition + case.setup.

    The standalone solver refuses npx*npy*npz > 1, and the cases spell their
    decomposition three different ways: three separate `par.nx = 2` lines
    (tpv8, tpv10, tpv104), one tuple `par.nx, par.ny, par.nz = 2, 2, 1`
    (drv.a6, tpv29), or not at all (tpv1053d, which inherits
    defaultParameters). APPENDING the override handles all three by
    construction -- it is the last assignment executed either way -- which
    matters because reading the decomposition instead of overriding it has
    twice left a case multi-rank and made the solver's refusal look like a
    coverage gap. par.HPC_ncpu keeps its original value; it only feeds the
    generated run.sh, which the standalone never reads.

    No Fortran binary is involved: case.setup writes the bFile/netCDF case
    inputs the standalone reads natively.
    """
    rc = subprocess.call([sys.executable,
                          os.path.join(REPO_ROOT, 'scripts', 'create.newcase'),
                          case_dir, case_name], env=env)
    if rc != 0:
        raise RuntimeError('create.newcase %s exited %d' % (case_name, rc))
    params = os.path.join(case_dir, 'user_defined_params.py')
    with open(params) as f:
        text = f.read().rstrip('\n')
    with open(params, 'w') as f:
        f.write(text + '\n\n# forced serial by testsys/e2e/run_e2e.py'
                       ' (the standalone solver is serial-only)\n'
                       'par.nx = 1\npar.ny = 1\npar.nz = 1\n')
    rc = subprocess.call([sys.executable, 'case.setup'], cwd=case_dir, env=env)
    if rc != 0:
        raise RuntimeError('case.setup for %s exited %d' % (case_name, rc))


def run_standalone(case_dir, backend, device='cpu', env=None):
    """`python3 -m eqdyna <case_dir> --backend <numpy|jax>` -- the
    exact command a user would type, with the backend ALWAYS named.

    It used to be invoked with no --backend at all, which silently meant the
    jax default: ~1700 lines of NumPy port (port.py, port_rsf.py, port_tp.py)
    had never been exercised by any gate while the tier reported success.
    """
    solver = {'python-numpy': 'numpy', 'python-jax': 'jax'}[backend]
    env = dict(env or base_env())
    env['PYTHONPATH'] = os.path.join(REPO_ROOT, 'src', 'python')
    env['PYTHONUNBUFFERED'] = '1'
    # JAX_PLATFORMS pins the device: a run labelled jax-on-cpu that silently
    # landed on a contended GPU is a different measurement under the same name.
    env['JAX_PLATFORMS'] = device
    rc = subprocess.call([sys.executable, '-u', '-m', 'eqdyna',
                          case_dir, '--backend', solver],
                         cwd=REPO_ROOT, env=env)
    if rc != 0:
        raise RuntimeError('eqdyna.eqdyna3d %s --backend %s exited %d'
                           % (os.path.basename(case_dir), solver, rc))
    frt = os.path.join(case_dir, 'frt.txt0')
    if not os.path.isfile(frt):
        raise RuntimeError('%s was not written by eqdyna.eqdyna3d' % frt)
    return frt


def run_fortran(case_name, case_dir, eqdyna_cmd, env):
    """create.newcase -> case.setup -> mpirun -> plotRuptureDynamics.

    plotRuptureDynamics is what writes fault.dyna.r.nc, the second artifact
    this backend is compared on (matrix.ARTIFACTS); dropping it once already
    turned that comparison into a comparison of a stale file.
    """
    test_dir = os.path.dirname(case_dir)
    steps = [
        (['create.newcase', os.path.basename(case_dir), case_name], test_dir),
        (['./case.setup'], case_dir),
        ([MPIRUN, '-np', str(matrix.FORTRAN_RANKS[case_name]), eqdyna_cmd], case_dir),
        ([sys.executable, 'plotRuptureDynamics'], case_dir),
    ]
    for cmd, cwd in steps:
        rc = _run(cmd, cwd, env)
        if rc != 0:
            raise RuntimeError('`%s` exited %d' % (' '.join(cmd), rc))


def run_cell(case, backend, test_dir, eqdyna_cmd, env, device):
    """Run one cell and return its run directory. Raises on any failure."""
    if backend == 'fortran':
        case_dir = os.path.join(test_dir, case)
        run_fortran(case, case_dir, eqdyna_cmd, env)
        return case_dir
    case_dir = os.path.join(test_dir, '%s.%s' % (case, backend))
    make_serial_case(case, case_dir, env)
    run_standalone(case_dir, backend, device=device, env=env)
    return case_dir


# --------------------------------------------------------------------------
# selection
# --------------------------------------------------------------------------
def select(args):
    """(runnable, declared_unsupported, label, explicit) for this invocation."""
    if args.ci:
        # --backends/--cases, WHEN COMBINED WITH --ci, filter matrix.CI_CELLS
        # itself rather than switching to matrix.cells() -- this is what lets
        # a CI matrix job ask for "the fortran slice of CI_CELLS" or "the
        # python slice of CI_CELLS" while matrix.CI_CELLS stays the ONE
        # declared, memory-measured source of truth (rule 6). The workflow
        # never spells out a cell list of its own; it only names an axis
        # subset, so widening CI_CELLS in matrix.py widens every job that
        # asks for 'all' on that axis with no workflow edit needed.
        wanted = list(matrix.CI_CELLS)
        if args.cases:
            want_cases = set(args.cases.split(','))
            wanted = [c for c in wanted if c[0] in want_cases]
        if args.backends:
            want_backends = set(args.backends.split(','))
            wanted = [c for c in wanted if c[1] in want_backends]
        runnable, unsupported = [], []
        for cell in wanted:
            if matrix.is_supported(*cell):
                runnable.append(cell)
            else:
                unsupported.append(cell + (matrix.unsupported_reason(*cell),))
        filt = ''
        if args.cases or args.backends:
            filt = (', filtered to cases=%s backends=%s'
                    % (args.cases or 'all', args.backends or 'all'))
        return (runnable, unsupported,
                'CI (declared cell list, chosen against a measured %.0f GB '
                'runner -- matrix.CI_CELLS%s)' % (matrix.CI_RUNNER_RAM_GB, filt),
                True)
    cases = args.cases.split(',') if args.cases else None
    backends = args.backends.split(',') if args.backends else None
    runnable, unsupported = matrix.cells(cases, backends)
    explicit = bool(cases or backends)
    label = ('explicit: cases=%s backends=%s'
             % (args.cases or 'all', args.backends or 'all')) if explicit \
        else 'default: every cell of the table'
    return runnable, unsupported, label, explicit


def memory_note(runnable):
    lines = ['measured peak RSS per cell (rule 6 -- the number travels with '
             'the decision):']
    for cell in runnable:
        rss = matrix.MEASURED_PEAK_RSS_GB.get(cell)
        lines.append('  %-16s %-13s %s'
                     % (cell[0], cell[1],
                        ('%.2f GB' % rss) if rss is not None
                        else 'not measured'))
    return lines


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument('--cases', help='comma-separated subset of the case axis')
    ap.add_argument('--backends', help='comma-separated subset of the backend axis')
    ap.add_argument('--ci', action='store_true',
                    help='run matrix.CI_CELLS, the declared memory-measured '
                         'cell list a 7 GB GitHub runner can hold')
    ap.add_argument('--jobs', type=int, default=None,
                    help='core budget for concurrent cells (default: cores-4). '
                         'A fortran cell costs its rank count, a python cell 1. '
                         '--jobs 1 runs serially, which is what a 2-core CI '
                         'runner should use; the printed output is identical '
                         'either way.')
    ap.add_argument('--device', default='cpu', choices=('cpu', 'cuda'),
                    help='JAX_PLATFORMS for the python-jax backend (default cpu)')
    args = ap.parse_args(argv)

    runnable, unsupported, label, explicit = select(args)

    print('\n==== e2e sweep: coverage ====')
    for line in matrix.coverage_report(runnable, unsupported, label):
        print(line)
    for line in memory_note(runnable):
        print(line)
    if args.device != 'cpu':
        print('device: JAX_PLATFORMS=%s for python-jax cells' % args.device)

    if explicit and unsupported:
        print('\ne2e: FAIL - the selection names %d cell(s) this table '
              'declares unsupported. A cell you asked for and could not be '
              'run is a failure, not a skip (rule 2):' % len(unsupported))
        for c, b, reason in unsupported:
            print('  - %s x %s: %s' % (c, b, reason))
        return 1
    if not runnable:
        print('\ne2e: FAIL - the selection runs zero cells; a sweep that '
              'compares nothing must never exit green')
        return 1

    # Gate 1 - cheap check before the expensive runs (rule 9).
    guard = os.path.join(REPO_ROOT, 'testsys', 'regression', 'test_create_newcase.py')
    if subprocess.call([sys.executable, guard]) != 0:
        print('e2e: FAIL - create.newcase guard failed; aborting before expensive runs')
        return 1

    # Gate 2 - fresh build, if and only if a fortran cell is selected.
    eqdyna_cmd = None
    if any(b == 'fortran' for _, b in runnable):
        if BIN_OVERRIDE:
            eqdyna_cmd = os.path.abspath(BIN_OVERRIDE)
            if not os.path.exists(eqdyna_cmd):
                print('e2e: FAIL - EQDYNA_E2E_BIN=%s does not exist' % BIN_OVERRIDE)
                return 1
            print('e2e: EQDYNA_E2E_BIN set - NOT rebuilding; this run says '
                  'nothing about the build-and-install path (rule 16). Using %s'
                  % eqdyna_cmd)
        else:
            bin_exe = os.path.join(REPO_ROOT, 'bin', 'eqdyna')
            if os.path.exists(bin_exe):
                os.remove(bin_exe)
            rc = subprocess.call(['./install-eqdyna.sh', '-m', MACHINE], cwd=REPO_ROOT)
            if rc != 0:
                print('e2e: FAIL - ./install-eqdyna.sh -m %s exited %d' % (MACHINE, rc))
                return 1
            if not os.path.exists(bin_exe):
                print('e2e: FAIL - build finished but bin/eqdyna does not exist')
                return 1
            eqdyna_cmd = 'eqdyna'
    else:
        print('e2e: no fortran cell in this selection - no Fortran build needed')

    # Rule 8 - preserve, never delete, the previous run's evidence.
    test_dir = os.path.join(REPO_ROOT, 'test')
    prev_dir = os.path.join(REPO_ROOT, 'test.prev')
    if os.path.isdir(test_dir):
        if os.path.isdir(prev_dir):
            shutil.rmtree(prev_dir)
        shutil.move(test_dir, prev_dir)
        print('e2e: preserved previous run as %s' % prev_dir)
    os.makedirs(test_dir)

    env = base_env()
    start = time.time()

    # Cells are independent: each gets its own directory ('<case>.<backend>'),
    # reads the same read-only reference tree, and shares nothing else. So the
    # sweep runs them concurrently rather than one at a time.
    #
    # WHY THIS MATTERS: measured serially, the 20-cell sweep took 3496 s on a
    # 64-core box, and 67% of that was the python-numpy column alone
    # (2267 s over 6 cells; test.tpv29 964 s, test.drv.a6 553 s). Meanwhile at
    # most 4 cores were busy. Wall time is now bounded by the SLOWEST CELL, not
    # the sum.
    #
    # The budget is CORES, not memory: all 20 cells together peak at ~23.5 GB
    # (matrix.MEASURED_PEAK_RSS_GB) against 1 TB here. A fortran cell costs its
    # rank count (4); a python cell is one process. --jobs sets the budget and
    # defaults to the machine's cores less a small reserve; --jobs 1 restores
    # the serial order exactly, which is what CI uses when its runner has 2.
    #
    # Output is COLLECTED per cell and printed when that cell finishes, never
    # streamed, so concurrent cells cannot interleave their lines into an
    # unreadable log. The results table is re-sorted into table order
    # afterwards, so a parallel run and a serial run print identically.
    def cell_cost(case, backend):
        if backend == 'fortran':
            # its real rank count from testNameList.coreNumList, not a guess
            return max(1, matrix.FORTRAN_RANKS.get(case, 4))
        return 1

    cells = [(c, b) for c in matrix.CASES for b in matrix.BACKENDS
             if (c, b) in runnable]
    order = {cb: i for i, cb in enumerate(cells)}

    def run_one(cb):
        case, backend = cb
        t0 = time.time()
        try:
            case_dir = run_cell(case, backend, test_dir, eqdyna_cmd, env,
                                args.device)
            ok, lines = compare.compare_cell(case, backend, case_dir)
        except Exception as exc:                # noqa: BLE001 - reported, not swallowed
            ok, lines = False, ['%s: %s' % (type(exc).__name__, exc)]
        return (case, backend, ok, time.time() - t0, lines)

    budget = args.jobs if args.jobs else max(1, (os.cpu_count() or 4) - 4)
    results = []
    if budget <= 1:
        for cb in cells:
            print('\n-- cell: %s x %s --' % cb)
            r = run_one(cb)
            results.append(r)
            for line in r[4]:
                print('   ' + line)
            print('%s %s x %s (%.1fs)'
                  % ('SUCCESS' if r[2] else 'FAIL', r[0], r[1], r[3]))
    else:
        import concurrent.futures as _cf
        import threading
        lock = threading.Lock()

        # ALL-OR-NOTHING core reservation. A cell costs more than one core
        # (a fortran cell costs its rank count), and the obvious spelling --
        # `for _ in range(cost): sem.acquire()` on a Semaphore -- DEADLOCKS,
        # because it takes units one at a time: with budget 6, two cells each
        # needing 4 can end up holding 3 apiece and both wait forever for a
        # fourth the other is holding. That is not hypothetical; it hung this
        # sweep for 80 minutes with 14 cells left, parent alive at 0.1% CPU
        # and no children. A Condition lets a thread take its whole cost
        # atomically or not at all, which cannot deadlock.
        cond = threading.Condition()
        free_cores = budget

        def reserve(cost):
            nonlocal free_cores
            with cond:
                while free_cores < cost:
                    cond.wait()
                free_cores -= cost

        def release(cost):
            nonlocal free_cores
            with cond:
                free_cores += cost
                cond.notify_all()

        def guarded(cb):
            # min(): a cell that costs more than the whole budget would
            # otherwise wait forever for cores that will never exist.
            cost = min(cell_cost(*cb), budget)
            reserve(cost)
            try:
                return run_one(cb)
            finally:
                release(cost)

        print('e2e: running %d cells concurrently, core budget %d '
              '(serial would be the sum of all cell times)'
              % (len(cells), budget))
        with _cf.ThreadPoolExecutor(max_workers=len(cells)) as pool:
            futs = {pool.submit(guarded, cb): cb for cb in cells}
            for fut in _cf.as_completed(futs):
                r = fut.result()
                with lock:
                    print('\n-- cell: %s x %s (%.1fs) --' % (r[0], r[1], r[3]))
                    for line in r[4]:
                        print('   ' + line)
                    print('%s %s x %s'
                          % ('SUCCESS' if r[2] else 'FAIL', r[0], r[1]))
                    results.append(r)
        results.sort(key=lambda r: order[(r[0], r[1])])
    elapsed = time.time() - start

    print('\n==== e2e sweep: results ====')
    print('%-16s %-13s %-8s %8s  %s' % ('case', 'backend', 'verdict', 'seconds',
                                        'measure'))
    for case, backend, ok, dt, lines in results:
        print('%-16s %-13s %-8s %8.1f  %s'
              % (case, backend, 'SUCCESS' if ok else 'FAIL', dt,
                 lines[0] if lines else ''))
    failed = [(c, b) for c, b, ok, _, _ in results if not ok]

    print('\n==== e2e sweep: SUMMARY ====')
    print('selection : %s' % label)
    print('ran       : %d of %d cells in the %d case x %d backend table '
          '(%d passed, %d failed)'
          % (len(results), len(matrix.CASES) * len(matrix.BACKENDS),
             len(matrix.CASES), len(matrix.BACKENDS),
             len(results) - len(failed), len(failed)))
    print('not gated : %d declared-unsupported cell(s): %s'
          % (len(unsupported),
             ', '.join('%s x %s' % (c, b) for c, b, _ in unsupported) or 'none'))
    print('wall clock: %.1fs' % elapsed)
    if len(results) != len(runnable):
        print('e2e: FAIL - %d cell(s) were selected but %d produced a verdict; '
              'a cell that produced no verdict is a failure'
              % (len(runnable), len(results)))
        return 1
    if failed:
        print('e2e: FAIL - %d cell(s) failed: %s'
              % (len(failed), ', '.join('%s x %s' % cb for cb in failed)))
        return 1
    print('e2e: SUCCESS - %d/%d selected cells matched test.reference.results/ '
          '(this states its own scope; see the coverage block above)'
          % (len(results), len(runnable)))
    return 0


if __name__ == '__main__':
    if '--full' in sys.argv:
        # The full tier (SCEC spec resolution/duration, 16 ranks, report-only)
        # is a different gate shape, not a cell of this sweep -- see
        # run_e2e_full.py's module docstring.
        sys.path.insert(0, E2E_DIR)
        from run_e2e_full import main as main_full
        sys.exit(main_full())
    sys.exit(main())

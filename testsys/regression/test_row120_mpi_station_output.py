#! /usr/bin/env python3
"""
Regression guard (board row 120): python-jax-mpi's per-rank station
SELECTION (which faultst*/body* files each rank matches, via
eqdyna3d.build_solver_state -> meshgen.build_station_matching against THIS
RANK's local grid lines) reproduces Fortran's own per-rank selection
(setSurfaceStation/createMasterNode, meshgen.f90) EXACTLY, at test.tpv8's
opted-in decomposition -- including Fortran's two pre-existing gaps (rule
23: Fortran is the reference for numerics, warts included):
  - a station whose y sits exactly on a shared rank-boundary plane is
    dropped by every rank (no iy-edge branch anywhere in setSurfaceStation);
  - a station/fault-node on a shared x-boundary plane (which DOES have an
    edge branch) can be matched by more than one rank and therefore written
    under the same filename by more than one rank -- the same "last rank to
    call wins" shape testsys/matrix.py's GATE_STATIONS docstring already
    documents for Fortran's own multi-rank runs.

THE BEHAVIOUR under test (rule 10a: on behaviour, not a source grep):
  1. A fresh 4-rank Fortran run of test.tpv8 (tiny term, 5 steps -- this
     guard is about SELECTION, not physics) writes some SET of faultst*/
     body* files. That set is the ground truth.
  2. Calling eqdyna3d.build_solver_state(case_dir, part=Partition(rank, 4,
     (2, 2, 1))) for rank in 0..3 -- IN PROCESS, no mpirun, no jax, no
     time-stepping -- and taking the UNION, across ranks, of the filenames
     each rank's S['st_on_*']/S['st_off_*'] would produce (via
     library_output.onfault_filename/offfault_filename) reproduces THAT
     EXACT SET.
  3. Station 11 ([0, 0.5, -0.3] km, y=+500 m -- exactly tpv8's mey=0/mey=1
     boundary at this decomposition, docs/notes/NOTES_row120.md) is absent
     from BOTH sets: the known drop is reproduced, not silently "fixed" by
     one side and not by the other.
  4. eqdyna3d.run_case_mpi's ACTUAL writer (library_output.
     write_onfault_stations/write_offfault_stations, fed one rank's real
     S plus a synthetic one-step history of the right shape) writes exactly
     the filenames step 2 predicted for that rank -- so the filename
     FORMULA and the WRITER are not silently out of sync with each other.

MUTATION CHECK (test_station_gate.py's idiom: mutate the compared DATA, not
the source under test): the set-equality check that step 2 relies on is
exercised once on the real (matching) data (must PASS) and once on a
deliberately corrupted copy -- one filename dropped, one spurious filename
added (must FAIL) -- so a reviewer can see the check has teeth, not just
that it happened to pass once.

Cost: builds ONLY src/fortran/eqdyna (never bin/eqdyna, matching
test_fault_mpi_boundary_arn.py -- a shared box may have another session
running bin/eqdyna concurrently), runs it for 5 steps at 4 ranks on the
REAL test.tpv8 mesh (the geometry, not a shrunk stand-in, is the point: the
x=0/y=+500 rank-boundary coincidences only exist at tpv8's own dx). The
python side never launches mpirun or jax at all -- build_solver_state's
per-rank matching and the writer are both pure, MPI-free functions.
"""
import glob
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FSRC = os.path.join(ROOT, 'src', 'fortran')
PYSRC = os.path.join(ROOT, 'src', 'python')
sys.path.insert(0, PYSRC)

CASE = 'test.tpv8'
RANKS = 4
DECOMP = (2, 2, 1)
MACHINE = os.environ.get('EQDYNA_TEST_MACHINE', 'ubuntu')
MPIRUN = os.environ.get('EQDYNA_MPIRUN', 'mpirun')

# The one requested off-fault station (1-indexed, case_input/test.tpv8's
# st_coor_off_fault) this mission's evidence identified as sitting exactly
# on the mey=0/mey=1 boundary at (2,2,1): [0, 0.5, -0.3] km.
DROPPED_STATION_FILE = 'body005st000dp003.txt'


def _env():
    env = dict(os.environ)
    env['EQDYNAROOT'] = ROOT
    env['PATH'] = os.pathsep.join([os.path.join(ROOT, 'bin'),
                                   os.path.join(ROOT, 'scripts'),
                                   env.get('PATH', '')])
    return env


def build_eqdyna():
    """Build ONLY src/fortran/eqdyna -- never touch bin/eqdyna (a shared box
    may have another session running it concurrently)."""
    env = dict(os.environ)
    env['MACHINE'] = MACHINE
    r = subprocess.run(['make'], cwd=FSRC, env=env, capture_output=True,
                       text=True, timeout=300)
    binpath = os.path.join(FSRC, 'eqdyna')
    if r.returncode != 0 or not os.path.exists(binpath):
        raise RuntimeError(
            'src/fortran build failed (exit %d); is mpif90 (MACHINE=%s) on '
            'PATH?\n%s' % (r.returncode, MACHINE, (r.stdout + r.stderr)[-3000:]))
    return binpath


def make_case(case_dir):
    """create.newcase + a tiny par.term override (5 steps -- this guard is
    about station SELECTION, not physics) + case.setup. Leaves the case's
    own par.nx/ny/nz at test.tpv8's committed (2, 2, 1) (case_input/
    test.tpv8/user_defined_params.py), matching MPI4NodalQuant.DECOMP[4] --
    both the fresh Fortran binary and (per build_solver_state's own check)
    the python port accept this same decomposition."""
    env = _env()
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'create.newcase'),
                        case_dir, CASE], env=env, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError('create.newcase failed: %s' % (r.stderr or '')[-1500:])
    params = os.path.join(case_dir, 'user_defined_params.py')
    with open(params) as f:
        text = f.read().rstrip('\n')
    with open(params, 'w') as f:
        f.write(text + '\n\n# tiny term override by '
                       'test_row120_mpi_station_output.py (SELECTION, not '
                       'physics, is under test)\n'
                       'par.term = 5.0 * par.dt\n')
    r = subprocess.run([sys.executable, 'case.setup'], cwd=case_dir,
                       env=env, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError('case.setup failed: %s' % (r.stderr or '')[-1500:])


def _station_files(case_dir):
    names = set()
    for pat in ('faultst*.txt', 'body*.txt'):
        names |= {os.path.basename(p) for p in glob.glob(os.path.join(case_dir, pat))}
    return names


def run_fortran(binpath, case_dir):
    r = subprocess.run([MPIRUN, '-np', str(RANKS), binpath], cwd=case_dir,
                       capture_output=True, text=True, timeout=180)
    if r.returncode != 0:
        raise RuntimeError('fortran 4-rank run exited %d:\n%s'
                           % (r.returncode, (r.stdout + r.stderr)[-3000:]))
    return _station_files(case_dir)


def python_selection(case_dir):
    """Per-rank filename UNION, pure Python, no mpirun/jax/time-stepping --
    (rank_files, union) where rank_files[r] is rank r's own filename set
    (needed by the writer-agreement check) and union is their union."""
    from eqdyna import MPI4NodalQuant as MQ
    from eqdyna import eqdyna3d, library_output

    rank_files = []
    rank_S = []
    for rank in range(RANKS):
        part = MQ.Partition(rank, RANKS, DECOMP)
        S, _mesh = eqdyna3d.build_solver_state(case_dir, part=part)
        files = set()
        for i in range(int(S['st_on_idx'].shape[0])):
            files.add(library_output.onfault_filename(
                float(S['st_on_strike_m'][i]), float(S['st_on_depth_m'][i]),
                S['fault_dip_rad']))
        for i in range(int(S['st_off_idx'].shape[0])):
            files.add(library_output.offfault_filename(
                float(S['st_off_x_m'][i]), float(S['st_off_y_m'][i]),
                float(S['st_off_z_m'][i])))
        rank_files.append(files)
        rank_S.append(S)
    union = set()
    for f in rank_files:
        union |= f
    return rank_S, rank_files, union


def check_writer_agreement(rank_S, rank_files, scratch):
    """Feeds ONE rank's real S (whichever has both on- and off-fault
    matches, if any does) plus a synthetic one-step zero history into the
    ACTUAL writer eqdyna3d.run_case_mpi calls, and asserts the files that
    land on disk equal what python_selection predicted for that rank --
    catching a formula/writer mismatch that a pure-filename check could not."""
    import numpy as np
    from eqdyna import library_output

    for rank, S in enumerate(rank_S):
        n_on = int(S['st_on_idx'].shape[0])
        n_off = int(S['st_off_idx'].shape[0])
        if n_on == 0 and n_off == 0:
            continue
        rank_dir = os.path.join(scratch, 'rank%d' % rank)
        os.makedirs(rank_dir, exist_ok=True)
        ncol_on = 11 if S['friclaw'] >= 3 else 8
        on_hist = np.zeros((n_on, ncol_on, 1))
        off_hist = np.zeros((n_off, 7, 1))
        on_paths = library_output.write_onfault_stations(rank_dir, S, on_hist)
        off_paths = library_output.write_offfault_stations(rank_dir, S, off_hist)
        written = {os.path.basename(p) for p in on_paths + off_paths}
        if written != rank_files[rank]:
            raise AssertionError(
                'rank %d: writer produced %r, predicted (from the same S) %r'
                % (rank, sorted(written), sorted(rank_files[rank])))
        return True  # at least one rank checked
    return False


def _selection_matches(fortran_set, python_set):
    """The ONE comparison this guard is built on -- factored out so the
    MUTATION check below can call it on corrupted data too."""
    return fortran_set == python_set


def main():
    if shutil.which('mpif90') is None or shutil.which(MPIRUN) is None:
        print('FAIL test_row120_mpi_station_output: mpif90/%s not on PATH '
              '(this guard requires a real MPI build, it does not skip '
              'silently)' % MPIRUN)
        return 1

    fails = []
    try:
        binpath = build_eqdyna()
    except RuntimeError as e:
        print('FAIL test_row120_mpi_station_output')
        print(' -', e)
        return 1

    tmp = tempfile.mkdtemp(prefix='testRow120Station.')
    try:
        case_dir = os.path.join(tmp, CASE)
        make_case(case_dir)
        fortran_set = run_fortran(binpath, case_dir)
        rank_S, rank_files, python_set = python_selection(case_dir)

        checks = []

        checks.append(('fortran wrote at least one station file',
                       len(fortran_set) > 0, True))
        checks.append(('unmutated: fortran set == python union',
                       _selection_matches(fortran_set, python_set), True))
        checks.append((
            'MUTATION: python union missing one real file',
            _selection_matches(fortran_set, python_set - {next(iter(python_set))})
            if python_set else True,
            False))
        checks.append((
            'MUTATION: python union has one spurious file',
            _selection_matches(fortran_set, python_set | {'body999st999dp999.txt'}),
            False))
        checks.append((
            'station 11 ([0, 0.5, -0.3] km, y-partition boundary) is '
            'dropped by BOTH fortran and python',
            DROPPED_STATION_FILE not in fortran_set
            and DROPPED_STATION_FILE not in python_set,
            True))

        writer_checked = check_writer_agreement(rank_S, rank_files, tmp)
        checks.append(('writer output matches the predicted filenames for '
                       'at least one rank', writer_checked, True))

        for label, got, expect_pass in checks:
            ok = (got is True) == expect_pass
            print(('PASS' if ok else 'FAIL') + '  ' + label)
            if not ok:
                fails.append(label)

        print('fortran station files (%d): %s' % (len(fortran_set), sorted(fortran_set)))
        print('python  station files (%d): %s' % (len(python_set), sorted(python_set)))
    except (RuntimeError, AssertionError) as e:
        print('FAIL test_row120_mpi_station_output')
        print(' -', e)
        return 1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    n_mut = sum(1 for label, *_ in checks if label.startswith('MUTATION'))
    print('%s test_row120_mpi_station_output: %d check(s) (%d mutation(s)); '
          '%d failure(s)' % ('FAIL' if fails else 'SUCCESS', len(checks),
                              n_mut, len(fails)))
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())

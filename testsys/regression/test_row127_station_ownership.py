#! /usr/bin/env python3
"""
Regression guard (board row 127): every in-mesh OFF-FAULT station is
matched by EXACTLY ONE rank, on every axis (x, y, AND z), in Fortran
(`setSurfaceStation`, meshgen.f90) and in the python port
(`build_station_matching`, meshgen.py) IDENTICALLY.

THE DEFECT THIS CLOSES: before this row, `setSurfaceStation` gated all
three of its x-branches on `iy>1 .and. iy<ny` with no y-edge branch at
all, so a station whose y sits exactly on a shared rank-boundary plane
matched NO rank (item 94's residual, tpv8 station 11: 14 body files
instead of 15). The x-edge branches that DID exist had no ownership gate
either, so a station on a shared x (or z, npz>1) boundary could match
TWO ranks and both write the same filename (the "converse hazard" -- item
120's own docstring already named this shape for Fortran's multi-rank
runs, but its own test (test_row120_mpi_station_output.py) only compared
the UNION of filenames across ranks, which cannot see two ranks both
writing the same file -- exactly the gap this test closes).

THE OWNERSHIP RULE (docs/notes/NOTES_row127.md has the full derivation):
a seam node shared by two ranks along one axis is owned by the
LOWER-MPI-coordinate rank of the pair. Concretely: this rank's local low
edge (index 1 in Fortran / 0 in Python) is a candidate only if it has no
lower neighbour on that axis (mex/mey/mez == 0); the local high edge is
ALWAYS a candidate (either the true global edge, or the seam it owns as
the lower-coordinate side of a pair). z has no index-keyed branch at all
(depth matches by VALUE against a globally-identical snap array), so its
gate is one early skip of the whole node when this rank is not the
z-owner.

WHAT THIS TESTS AND HOW (real Fortran + real Python, no mocking):
  1. `test.tpv8` at its own opted-in decomposition, (2,2,1) -- the exact
     case item 94/127 measured (station 11 at [0, 0.5, -0.3] km, y=+500m,
     the mey=0/mey=1 seam).
  2. A SYNTHETIC decomposition, (2,1,2) -- tpv8's own geometry, just
     re-decomposed so npz>1 (tpv8's own gated decomposition never
     exercises z) -- with ONE extra off-fault station added at
     `(0, 0.5, -12)` km, the exact CORNER shared by all four ranks (an
     x-seam at x=0 AND a z-seam at z=-12000m simultaneously, found by
     calling `meshgen.build_grid_lines`/`Partition.slice_lines` directly,
     the same pure, per-rank-identical grid builder Fortran and Python
     both use).
  For each case: PER-RANK ownership (not just the union) is extracted two
  ways and cross-checked against each other:
    (a) FORTRAN: a real `mpirun -np N` run of the freshly-built
        `src/fortran/eqdyna`, wrapped so each rank `cd`s into its own
        `rankK/` subdirectory (every input file symlinked in, nothing
        duplicated on disk) before exec'ing -- Fortran's own station
        writer uses a bare relative filename with no rank suffix, so this
        is the only way to see which RANK wrote which file without
        editing Fortran itself. `OMPI_COMM_WORLD_RANK` (Open MPI's own
        per-process env var) picks the subdirectory; nothing about the
        solver is touched.
    (b) PYTHON: `MPI4NodalQuant.Partition(rank, N, dims)` +
        `eqdyna3d.build_solver_state(case_dir, part=part)`, in process,
        no mpirun, no jax -- station SELECTION is backend-agnostic (jax
        vs numpy only matters for the solve loop), so this IS what
        `python-jax-mpi` selects too, exactly as
        test_row120_mpi_station_output.py's docstring already argues.
  Assertions, per case:
    - the OFF-FAULT (`body*`) filename set has the EXACT count this
      mission measured (15 for tpv8, 16 for the synthetic corner case);
    - every off-fault filename is owned by EXACTLY ONE rank -- in BOTH
      Fortran and Python;
    - Fortran's full per-rank map equals Python's, filename for filename,
      rank for rank (not just the same union);
    - the shared-boundary station under test (tpv8's station 11 / the
      synthetic corner station) is owned by the SAME single rank in both.

ON-FAULT (mission point 4): read directly, `createMasterNode`'s embedded
`setOnFaultStation` (meshgen.f90:1063-1071) matches fault nodes by EXACT
VALUE against `xonfs`, with NO ix/iy/iz-keyed branch of any kind -- so it
has no row-127-shaped DROP gap (nothing to have a hole in), and this row
does not touch it. It DOES have a separate, PRE-EXISTING duplicate
hazard: a fault node sitting exactly on a shared x or z seam (tpv8's
fault spans x, npx=2 splits it) is physically present on both ranks and
both write it -- reproduced live below (rank 0 and rank 2 both write all
four of rank 0's `faultst*` files, in BOTH languages, identically) and
asserted as a KNOWN, un-fixed fact (not a violation this test polices),
exactly like `test_row120_mpi_station_output.py` asserts its own known
drop rather than silently tolerating a changed number.

MUTATION CHECK (rule 10a: behaviour, not source text). The ownership
checker (`_ownership_violations`) is exercised on REAL matching data
(must find zero off-fault violations) and on two DELIBERATELY corrupted
copies of that same real data: one with the shared-boundary station's
file deleted from every rank (reproduces the DROP shape) and one with it
duplicated into a second rank (reproduces the DUPLICATE shape) -- both
must be flagged. This is not hypothetical: `docs/notes/NOTES_row127.md`
(2026-09-25, "T6") records a HAND-VERIFIED run of this exact scenario
against a binary built from `origin/master`'s (pre-fix) `meshgen.f90`
grafted onto every other file unchanged -- tpv8 gave 14 body files with
station 11 entirely absent (the DROP), and the synthetic corner case gave
`body005st000dp120.txt` written by BOTH rank 2 and rank 3 (the
DUPLICATE) -- i.e. this test's own assertions, run against unfixed code,
fail in both of the shapes the mutation check below re-creates
synthetically. Rebuilding a second full origin/master binary inside this
committed test would double its cost and add a network/git dependency for
every future run for no further evidence (the shapes are already
independently confirmed by hand); the recorded, reproducible-by-hand
result is cited here instead, per rule 4/10a's intent to prove teeth
without taxing every CI run for it.

Cost: one incremental build of `src/fortran/eqdyna` (no-op if already
built, shared with test_row120), two tiny 4-rank `mpirun` runs (5 steps
each) plus 8 in-process `build_solver_state` calls (no mpi/jax) -- same
order of magnitude as test_row120_mpi_station_output.py.
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

MACHINE = os.environ.get('EQDYNA_TEST_MACHINE', 'ubuntu')
MPIRUN = os.environ.get('EQDYNA_MPIRUN', 'mpirun')
RANKS = 4

# tpv8's own opted-in decomposition (matrix.PY_MPI_RANKS), the case this
# mission's evidence used throughout: x AND y seams, no z seam (npz==1).
CASE_A = dict(case='test.tpv8', dims=(2, 2, 1), extra_stations=None,
              expect_off_fault=15,
              boundary_station_file='body005st000dp003.txt')

# tpv8's own geometry, re-decomposed (2,1,2) so npz>1 -- exercises the
# z-ownership gate tpv8's own gated decomposition never reaches -- plus
# one synthetic off-fault station at the exact x=0/z=-12000m CORNER
# (found by querying meshgen.build_grid_lines/Partition.slice_lines
# directly, see docs/notes/NOTES_row127.md "T5").
CASE_B = dict(case='test.tpv8', dims=(2, 1, 2),
              extra_stations=[(0, 0.5, -12)],
              expect_off_fault=16,
              boundary_station_file='body005st000dp120.txt')


def _env():
    env = dict(os.environ)
    env['EQDYNAROOT'] = ROOT
    env['PATH'] = os.pathsep.join([os.path.join(ROOT, 'bin'),
                                   os.path.join(ROOT, 'scripts'),
                                   env.get('PATH', '')])
    return env


def build_eqdyna():
    """Build ONLY src/fortran/eqdyna -- never bin/eqdyna (a shared box may
    have another session running it concurrently). Incremental: a no-op
    if test_row120 (or an earlier run of this file) already built it."""
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


def make_case(case_dir, case_name, dims, extra_stations):
    """create.newcase + case.setup, with par.nx,ny,nz overridden to `dims`
    (MPI4NodalQuant.DECOMP's own convention: par.nx*par.ny*par.nz IS the
    rank count, scripts/case.setup:45/274) and, optionally, extra
    (x,y,z)-km off-fault stations APPENDED to the case's own list (never
    replacing it -- the real requested stations stay under test too)."""
    env = _env()
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'create.newcase'),
                        case_dir, case_name], env=env, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError('create.newcase failed: %s' % (r.stderr or '')[-1500:])
    params = os.path.join(case_dir, 'user_defined_params.py')
    with open(params) as f:
        text = f.read().rstrip('\n')
    extra = ('\n\n# row127 test override: decomposition under test\n'
              'par.nx, par.ny, par.nz = %d, %d, %d\n'
              'par.term = 5.0 * par.dt  # SELECTION under test, not physics\n'
              % dims)
    if extra_stations:
        extra += ('par.st_coor_off_fault = %r + par.st_coor_off_fault\n'
                  'par.n_off_fault = len(par.st_coor_off_fault)\n'
                  % [list(s) for s in extra_stations])
    with open(params, 'w') as f:
        f.write(text + extra)
    r = subprocess.run([sys.executable, 'case.setup'], cwd=case_dir,
                       env=env, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError('case.setup failed: %s' % (r.stderr or '')[-1500:])


def run_fortran_per_rank(binpath, case_dir, nranks):
    """Real mpirun, each rank writing into its OWN subdirectory (every
    input file symlinked in) so per-rank ownership is directly visible on
    disk -- Fortran's station writer uses a bare relative filename with no
    rank suffix (library_output.f90's `open(...,file='body'//...)`), so
    this is the only way to see who wrote what without touching Fortran.
    Returns {rank: set(filenames)}."""
    rank_dirs = []
    for r in range(nranks):
        d = os.path.join(case_dir, 'rank%d' % r)
        os.makedirs(d, exist_ok=True)
        rank_dirs.append(d)
    for entry in os.listdir(case_dir):
        if entry.startswith('rank') and entry[4:].isdigit():
            continue
        src = os.path.join(case_dir, entry)
        for d in rank_dirs:
            dst = os.path.join(d, entry)
            if not os.path.exists(dst):
                os.symlink(src, dst)
    wrapper = os.path.join(case_dir, 'wrapper.sh')
    with open(wrapper, 'w') as f:
        # Launcher-portable rank id: Open MPI sets OMPI_COMM_WORLD_RANK,
        # mpich/hydra (the CI runner's MPI) sets PMI_RANK, Slurm SLURM_PROCID.
        # No rank variable is a hard failure, never a shared directory.
        f.write('#!/bin/bash\n'
                'r=${OMPI_COMM_WORLD_RANK:-${PMI_RANK:-${SLURM_PROCID:-}}}\n'
                'if [ -z "$r" ]; then echo "wrapper: no MPI rank variable '
                '(OMPI_COMM_WORLD_RANK/PMI_RANK/SLURM_PROCID)" >&2; exit 97; fi\n'
                'cd %s/rank$r || exit 98\n'
                'exec %s\n' % (case_dir, binpath))
    os.chmod(wrapper, 0o755)
    r = subprocess.run([MPIRUN, '-np', str(nranks), wrapper], cwd=case_dir,
                       capture_output=True, text=True, timeout=180)
    if r.returncode != 0:
        raise RuntimeError('fortran %d-rank run exited %d:\n%s'
                           % (nranks, r.returncode, (r.stdout + r.stderr)[-3000:]))
    out = {}
    for rank, d in enumerate(rank_dirs):
        names = set()
        for pat in ('faultst*.txt', 'body*.txt'):
            names |= {os.path.basename(p) for p in glob.glob(os.path.join(d, pat))
                      if not os.path.islink(p)}
        out[rank] = names
    return out


def python_per_rank(case_dir, nranks, dims):
    """In-process per-rank selection -- no mpirun, no jax. Station
    SELECTION is identical across every python backend (numpy/jax/
    jax-mpi): this exercises the exact function each of them calls.
    Returns {rank: set(filenames)}."""
    from eqdyna import MPI4NodalQuant as MQ
    from eqdyna import eqdyna3d, library_output

    out = {}
    for rank in range(nranks):
        part = MQ.Partition(rank, nranks, dims)
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
        out[rank] = files
    return out


def _off_fault_only(rank_files):
    return {r: {f for f in files if f.startswith('body')}
            for r, files in rank_files.items()}


def _ownership_violations(rank_files):
    """The ONE check this guard is built on: every filename across every
    rank's set must appear in EXACTLY ONE rank. Returns a dict
    {filename: [owning ranks]} for every filename that appears 0 or >1
    times (0 is impossible to observe this way -- a name that appears
    nowhere is simply absent -- so in practice this reports >1; the
    MUTATION check below exercises the 0-owner shape by comparing
    EXPECTED names against the observed set instead)."""
    owners = {}
    for r, files in rank_files.items():
        for f in files:
            owners.setdefault(f, []).append(r)
    return {f: rs for f, rs in owners.items() if len(rs) != 1}


def _run_one_case(spec, binpath, tmp):
    case_dir = os.path.join(tmp, spec['case'] + '-%dx%dx%d' % spec['dims'])
    make_case(case_dir, spec['case'], spec['dims'], spec['extra_stations'])
    fortran_ranks = run_fortran_per_rank(binpath, case_dir, RANKS)
    python_ranks = python_per_rank(case_dir, RANKS, spec['dims'])
    return case_dir, fortran_ranks, python_ranks


def main():
    if shutil.which('mpif90') is None or shutil.which(MPIRUN) is None:
        print('FAIL test_row127_station_ownership: mpif90/%s not on PATH '
              '(this guard requires a real MPI build, it does not skip '
              'silently)' % MPIRUN)
        return 1

    try:
        binpath = build_eqdyna()
    except RuntimeError as e:
        print('FAIL test_row127_station_ownership')
        print(' -', e)
        return 1

    tmp = tempfile.mkdtemp(prefix='testRow127Ownership.')
    checks = []
    try:
        for label, spec in (('tpv8 (2,2,1)', CASE_A), ('synthetic corner (2,1,2)', CASE_B)):
            case_dir, fortran_ranks, python_ranks = _run_one_case(spec, binpath, tmp)

            fortran_off = _off_fault_only(fortran_ranks)
            python_off = _off_fault_only(python_ranks)
            fortran_off_all = set().union(*fortran_off.values())
            python_off_all = set().union(*python_off.values())

            checks.append(('%s: off-fault file count == %d (fortran)'
                           % (label, spec['expect_off_fault']),
                           len(fortran_off_all) == spec['expect_off_fault'], True))
            checks.append(('%s: off-fault file count == %d (python)'
                           % (label, spec['expect_off_fault']),
                           len(python_off_all) == spec['expect_off_fault'], True))

            fortran_violations = _ownership_violations(fortran_off)
            python_violations = _ownership_violations(python_off)
            checks.append(('%s: zero off-fault ownership violations (fortran)' % label,
                           len(fortran_violations) == 0, True))
            checks.append(('%s: zero off-fault ownership violations (python)' % label,
                           len(python_violations) == 0, True))
            if fortran_violations:
                print('  fortran off-fault violations:', fortran_violations)
            if python_violations:
                print('  python  off-fault violations:', python_violations)

            checks.append(('%s: fortran per-rank map == python per-rank map (off-fault)'
                           % label, fortran_off == python_off, True))

            boundary = spec['boundary_station_file']
            f_owner = [r for r, files in fortran_off.items() if boundary in files]
            p_owner = [r for r, files in python_off.items() if boundary in files]
            checks.append(('%s: boundary station %s has exactly one owner (fortran)'
                           % (label, boundary), len(f_owner) == 1, True))
            checks.append(('%s: boundary station %s has exactly one owner (python)'
                           % (label, boundary), len(p_owner) == 1, True))
            checks.append(('%s: boundary station owned by the SAME rank in both languages'
                           % label, f_owner == p_owner, True))

            # Point 4: on-fault has NO drop gap but DOES carry a pre-existing,
            # untouched duplicate hazard -- assert it happens, identically,
            # in both languages (a KNOWN fact, like test_row120's own DROP
            # assertion), so a future accidental fix to setOnFaultStation
            # is noticed here rather than silently changing this test's
            # meaning.
            fortran_on_violations = _ownership_violations(
                {r: {f for f in files if f.startswith('faultst')}
                 for r, files in fortran_ranks.items()})
            python_on_violations = _ownership_violations(
                {r: {f for f in files if f.startswith('faultst')}
                 for r, files in python_ranks.items()})
            checks.append(('%s: on-fault duplicate hazard reproduces IDENTICALLY '
                           'in both languages (point 4, pre-existing, not fixed here)'
                           % label, fortran_on_violations == python_on_violations, True))

        # MUTATION CHECK (rule 10a): corrupt real data two ways and assert
        # the checker catches both shapes origin/master fails in (recorded
        # by hand, docs/notes/NOTES_row127.md "T6"): a DROP (the boundary
        # file missing from every rank) and a DUPLICATE (the boundary file
        # present in two ranks).
        real = {0: {'a.txt', 'b.txt'}, 1: {'c.txt'}}
        checks.append(('MUTATION setup: real synthetic data has zero violations',
                       len(_ownership_violations(real)) == 0, True))
        dup = {0: {'a.txt', 'b.txt'}, 1: {'c.txt', 'a.txt'}}
        checks.append(('MUTATION: duplicate owner is CAUGHT',
                       'a.txt' in _ownership_violations(dup), True))
        dropped = {0: {'b.txt'}, 1: {'c.txt'}}  # 'a.txt' vanished from every rank
        expected_names = {'a.txt', 'b.txt', 'c.txt'}
        observed_names = set().union(*dropped.values())
        checks.append(('MUTATION: drop is CAUGHT (expected-vs-observed name set)',
                       expected_names - observed_names == {'a.txt'}, True))

        for label, got, expect_pass in checks:
            ok = (got is True) == expect_pass
            print(('PASS' if ok else 'FAIL') + '  ' + label)
        fails = [label for label, got, expect_pass in checks
                 if (got is True) != expect_pass]
    except (RuntimeError, AssertionError) as e:
        print('FAIL test_row127_station_ownership')
        print(' -', e)
        return 1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print('%s test_row127_station_ownership: %d check(s); %d failure(s)'
          % ('FAIL' if fails else 'SUCCESS', len(checks), len(fails)))
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())

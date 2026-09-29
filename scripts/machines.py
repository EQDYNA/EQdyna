#! /usr/bin/env python3
"""
The one HPC/machine registry (owner request 2026-09-29: run.py should sbatch
the test suite with a machine flag, "a refactored way to add this support",
"you shouldn't reinvent everything").

Adding a new HPC means adding ONE entry to MACHINES below, plus its own
`install-eqdyna.sh` branch (module loads / venv setup / compiler selection
stay there -- this module never loads a module file or builds a venv itself,
it only describes the machine for scripts/case.setup and testsys/run.py).

Fields per machine:
  scheduler       'slurm' or None (no batch scheduler -- run in place).
  partition       default `#SBATCH -p` value, or None if the machine has no
                   sane default (submit then requires --partition).
  cores_per_node  cores in one node, or None if not yet verified for this
                   machine (submit then refuses rather than guessing -- rule
                   2: no placeholder data).
  walltime        default `#SBATCH -t` value (hh:mm:ss), or None.
  mpirun          the launcher testsys' own `<launcher> -np N ...` invocation
                   (EQDYNA_MPIRUN) uses, or None if unverified (submit/
                   --machine then refuse rather than guess -- rule 2). NOT
                   necessarily the launcher a case's own batch.hpc uses --
                   TACC's `ibrun` takes `-n`/`-o`, not `-np` (no `-n` runs on
                   the WHOLE allocation), so it cannot serve run_e2e's
                   `-np N` calls; measured on LS6 2026-09-29 (idev node,
                   Intel MPI): plain `mpirun -np N` runs every cell
                   (fortran and python-jax-mpi) to completion.
  notes           one line, for humans reading this table.

This module also owns `write_slurm_header`, the SBATCH header shared by
every batch script this repo writes: scripts/case.setup's per-case
`batch.hpc` (formerly a private copy in case.setup, `_write_slurm_header`)
and testsys/run.py's `--submit` sweep job. One writer, so the two never
drift apart the way the three per-case writers did before case.setup
consolidated them (see case.setup's own comment on `_write_slurm_header`).
"""

MACHINES = {
    'ls6': dict(
        scheduler='slurm', partition='development', cores_per_node=128,
        walltime='02:00:00', mpirun='mpirun',
        notes='TACC Lonestar6; install-eqdyna.sh -e ls6 builds the venv '
              '(jax, mpi4py against Intel MPI) on $WORK. mpirun measured '
              '2026-09-29 (idev node, EQDYNA_MPIRUN unset -> Intel mpirun): '
              'every fortran and python-jax-mpi cell ran to completion '
              'under `mpirun -np N`. Pass --partition normal for a '
              'production-queue submission (development is the default '
              'here for a quick check).'),
    'grace': dict(
        scheduler='slurm', partition=None, cores_per_node=None,
        walltime='02:00:00', mpirun=None,
        notes='Texas A&M Grace; partition/cores-per-node/mpirun not yet '
              'verified here -- pass --partition explicitly, fill in '
              'cores_per_node, and confirm the launcher (see ls6\'s '
              'ibrun-vs-mpirun note) before relying on --submit or '
              '--machine grace'),
    'ubuntu': dict(
        scheduler=None, partition=None, cores_per_node=None,
        walltime=None, mpirun='mpirun',
        notes='local / CI, mpich -- no scheduler, runs in place'),
    'macos': dict(
        scheduler=None, partition=None, cores_per_node=None,
        walltime=None, mpirun='mpirun',
        notes='local, Homebrew mpich -- no scheduler, runs in place'),
}


def machine(name):
    """The registry entry for `name`, or raise naming the known machines."""
    if name not in MACHINES:
        raise ValueError(
            "unknown machine %r (known: %s) -- add it to "
            "scripts/machines.py's MACHINES table and to install-eqdyna.sh"
            % (name, ', '.join(sorted(MACHINES))))
    return MACHINES[name]


def write_slurm_header(f, *, jobname, nnode, ncpu, queue, walltime, account,
                        email, output='a.eqdyna.log%j'):
    """Write the `#SBATCH` header common to every batch script this repo
    generates. Kept as one function so case.setup's per-case batch.hpc and
    testsys/run.py's --submit sweep job cannot drift apart the way the three
    in-case writers did before case.setup itself consolidated them."""
    f.write("#! /bin/bash" + "\n")
    f.write("#SBATCH -J " + str(jobname) + "\n")
    f.write("#SBATCH -o " + str(output) + "\n")
    f.write("#SBATCH -N " + str(nnode) + "\n")
    f.write("#SBATCH -n " + str(ncpu) + "\n")
    f.write("#SBATCH -p " + str(queue) + "\n")
    f.write("#SBATCH -t " + str(walltime) + "\n")
    f.write("#SBATCH -A " + str(account) + "\n")
    # An empty email means "no one to notify" (case.setup's own default,
    # defaultParameters.py's HPC_email = ""), not "notify the empty string" --
    # emitting `--mail-user=` with nothing after the `=` is a directive with
    # no recipient, and `--mail-type` with no `--mail-user` at all is equally
    # meaningless, so both are skipped together.
    if email:
        f.write("#SBATCH --mail-user=" + str(email) + "\n")
        f.write("#SBATCH --mail-type=begin" + "\n")
        f.write("#SBATCH --mail-type=end" + "\n")

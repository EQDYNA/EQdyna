# Running a Case

## Create a case

```
create.newcase $caseDirectoryName $predefinedCompset
cd $caseDirectoryName
```

`$predefinedCompset` is one of the pre-defined cases on the Benchmarks page
(e.g. `test.tpv8`), or your own copy of one with `user_defined_params.py`
edited to match your fault and loading -- see the Parameters page for every
setting that file can override.

`TPV<number>` follows the naming convention of the SCEC/USGS Spontaneous
Rupture Code Verification Project; see the Benchmarks page for what each
one tests.

## Configure and run

```
./case.setup
bash run.sh
```

`case.setup` reads `user_defined_params.py`, writes the `FE_*.txt` input
files the solver reads, and writes `run.sh` (and, on an HPC cluster with a
batch queue, a submittable job script) for the case's own settings.
`run.sh` clears any previous output, runs the Fortran solver under MPI, and
then produces `cRuptureDynamics.png`.

## On an HPC cluster (LS6)

`case.setup` also writes `batch.hpc`, a SLURM job script, and `./case.submit`
submits it:

```
source install-eqdyna.sh -c ls6        # modules; sbatch passes this environment on
create.newcase mytpv8 test.tpv8
cd mytpv8
# in user_defined_params.py: par.HPC_account = "<your allocation>"
./case.setup
./case.submit
```

`case.submit` refuses a job with no project account. Queue, wall time and
node count come from `HPC_queue`, `HPC_time` and `HPC_ncpu` (see
[Parameters](parameters.md)). The batch job runs the Fortran solver with
`mpirun -np <ranks>` (Intel MPI on LS6), then plots; per-case batch jobs do
not run the Python/JAX backend.

## MPI rank count

The number of MPI ranks is `nx * ny * nz` from `user_defined_params.py`
(the domain decomposition along each axis); `case.setup` writes that count
into `run.sh`'s `mpirun` line automatically, so changing `nx`/`ny`/`nz` and
re-running `./case.setup` is enough to change the rank count.

## Checking the mesh without solving

Both solvers screen the mesh right after it is built: each fault node must be
wired to cells on its own side of the fault, no cell may hold both nodes of a
split pair, fault-node areas must match the meshed facets, cell volumes and
Jacobians must be sane, and the fault-node frames must be orthonormal and
perpendicular to the fault. A mesh that fails stops the run with exit code
81 (sides), 82 (bridging), 83 (fault area), 84 (volume/mass) or 85 (normals),
and the `FATAL` block lists every check that failed; see
[Troubleshooting](troubleshooting.md).

`--mesh-only` runs the mesh build and this screen and then exits 0 without
solving or writing outputs. From the case directory, `mpirun -np 1 eqdyna
--mesh-only` runs the Fortran binary, `python3 -m eqdyna . --mesh-only` the
Python port (serial only). Either prints one `mesh-screen: PASS` or
`mesh-screen: FAIL` line with the measured quantities.

## The Python/JAX backend

`run.sh` always runs the Fortran solver. The same case can instead be run
through EQdyna's Python/JAX port (`pip install jax` first). The Python solver
runs on one process, so set the decomposition to one rank before
`case.setup`, then run it from inside the case directory with `PYTHONPATH`
set as in the install step:

```
printf '\npar.nx = 1\npar.ny = 1\npar.nz = 1\n' >> user_defined_params.py
./case.setup
python3 -m eqdyna . --backend jax
```

Add `--device cuda` to run on an NVIDIA GPU (after `pip install "jax[cuda12]"`).

JAX is the default backend, and its default device is the CPU. A NumPy
backend also exists and can be run by hand with `--backend numpy`, but only
the Fortran and JAX solvers are covered by the project's automatic
verification against the SCEC benchmarks.

# Running a Case

## Create a case

```
create.newcase $caseDirectoryName $predefinedCompset
cd $caseDirectoryName
```

`$predefinedCompset` is one of the pre-defined cases below, or your own
copy of one with `user_defined_params.py` edited to match your fault and
loading -- see the Parameters page for every setting that file can
override.

* test.drv.a6, for deterministic ground motion with a fractal fault and plasticity
* test.tpv8, test.tpv10, test.tpv36, test.tpv37, test.tpv104, test.tpv1053d, test.tpv29, test.tpv30

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

## MPI rank count

The number of MPI ranks is `nx * ny * nz` from `user_defined_params.py`
(the domain decomposition along each axis); `case.setup` writes that count
into `run.sh`'s `mpirun` line automatically, so changing `nx`/`ny`/`nz` and
re-running `./case.setup` is enough to change the rank count.

## The Python/JAX backend

`run.sh` always runs the Fortran solver. The same case can instead be run
through EQdyna's Python/JAX port, once `case.setup` has written its input
files, from inside the case directory:

```
python3 -m eqdyna . --backend jax
```

JAX is the default backend, and its default device is the CPU. A NumPy
backend also exists and can be run by hand with `--backend numpy`, but only
the Fortran and JAX solvers are covered by the project's automatic
verification against the SCEC benchmarks.

# LS6 (TACC Lonestar6) validation -- row 136

SLURM job `3480247`: `python3 testsys/run.py all --machine ls6 --submit` (the
`--submit` sbatch path, PR #52 `f9931a5`), run end-to-end for the first time on
a real SLURM system. Raw artifacts extracted at `<repo>/ls6/` from
`eqdyna_sweep_3480247.tgz`: `eqdyna_sweep_3480247.log` (42637 lines),
`.ledger.jsonl` (22 rows), `.profiles.jsonl` (55 rows).

## Environment

```
Currently Loaded Modules:
  1) intel/19.1.1    4) pmix/3.2.3   7) cmake/3.21.3
  2) impi/19.0.9     5) xalt/3.1     8) netcdf/4.6.2
  3) autotools/1.4   6) TACC         9) python/3.12.11
```
- `EQDYNAROOT=/scratch/07931/dunyuliu/eqdyna`
- eqdyna binary source stamp `333bcb5eb385` -- matches the tree at run time.
- Python venv: `/work/07931/dunyuliu/ls6/eqdyna-venv` (pytest 9.1.1).
- `testsys/run.py: --machine ls6 -> EQDYNA_TEST_MACHINE=ls6 EQDYNA_MPIRUN=mpirun`
  (Intel MPI's `mpirun`, not this dev box's Open MPI).

## Result

```
tiers run: unit, regression, e2e
SUCCESS unit (exit 0)
SUCCESS regression (exit 0)
FAIL e2e (exit 1)
run.py all --machine ls6 exit 1
```

- unit + regression: fully green.
- e2e: 22 of 23 gated cells pass. The one failure is `test.drv.a6 x fortran`
  (`FAIL nodes=5151 ruptured-in-both=1050 (matched-arrival=1014)`) -- this is
  the KNOWN ifort red on this chaotic case (row 138), which the owner ruled
  stays strict (gfortran-frozen references, frt flip gate passes 284/450).
  Not a new defect; expected going in.
- `test.drv.a6 x python-jax`: SUCCESS (189.6s) -- the jax backend passes this
  case on LS6 where fortran/ifort does not.
- Wall clock: 235.8s for the e2e sweep; 22 cell-wall-clock rows appended to
  `docs/perf_ledger.jsonl` (snapshot
  `docs/perf_snapshots/e2e_cells_2026-09-30_080311_2667377.json`).
- Station-series and nc-threshold gates passed at their normal bounds on every
  passing cell (e.g. `test.tpv36 x python-jax`: nc max|diff|=4.19e-09 vs bound
  1e-06; worst station e=9.46e-11 vs bound 1e-07).
- 10 cells declared UNSUPPORTED as designed (`python-jax-mpi` opt-in not
  extended to those cases) -- not gated, not a red.

## Conclusion for row 136

The `--submit` SLURM path is validated end-to-end on a real allocation: job
submission, module environment, Intel MPI `mpirun`, and the full unit +
regression + e2e pipeline all completed and produced the expected verdict
table, with the single failure being the already-known, owner-accepted
`drv.a6 x fortran` ifort red. Nothing here changes row 136's status beyond
confirming the `--submit` path itself works; the `drv.a6` question is
unchanged and remains owner-held (row 138).

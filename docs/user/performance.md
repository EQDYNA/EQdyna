# Performance

Every number below is measured and dated. Figures are refreshed
periodically; treat older ones as historical rather than a live guarantee
for your own hardware.

## Per-step solver cost

Single core, `test.tpv8`, fixed setup/compile time subtracted by difference
over two step counts. Measured 2026-09-16 with `python3 testsys/run.py perf`:

| engine | ms/step | vs fortran |
|---|---|---|
| fortran | 305.8 | 1.00 |
| jax-cpu | 191.3 | 0.63 (faster) |
| numpy | 529.5 | 1.73 (slower) |

The metric is per-step cost, not total wall clock, so that compiler
warm-up time cannot be mistaken for a change in the solver itself.

## Peak memory (RSS)

One case at a time, measured with `/usr/bin/time -v`:

| case | numpy | jax |
|---|---|---|
| test.tpv8 | 1.45 GB | 1.91 GB |
| test.tpv10 | -- | 3.23 GB |
| test.tpv104 | 2.29 GB | 3.42 GB |
| test.tpv1053d | -- | 3.95 GB |
| test.drv.a6 | -- | 9.57 GB |

## Full benchmark sweep wall clock

Between 2026-09-17 and 2026-09-19, the local sweep (`python3 testsys/run.py
all`) took between about 2100 and 2300 seconds across three runs on a 64-core
workstation. That sweep then ran 30 cells, including the NumPy backend and
each case at its own full simulated time; today's sweep runs fewer cells at a
5 s term and is faster, and has not been re-timed here. Wall clock depends heavily on
other work sharing the same machine at the time. Cases are independent and
each uses roughly one core, so throughput comes from running many cases at
once, not from any single case scaling across cores.

## Core scaling

Core-count scaling of a single JAX-CPU run is not currently characterised.
One ratio is measured on current code: 193 ms/step pinned to one core versus
13.8 ms/step unpinned, about 14x. An older curve (2.00x / 3.96x / 7.48x on
2 / 4 / 8 cores) predates a solver restructure and has not been reproduced.
Its knee sits at exactly 8 cores, the NUMA node size of the machine, so it may
measure memory locality rather than the solver. Treat any core-scaling claim
as provisional until it is remeasured.

## Large runs on an HPC cluster

Measured on Lonestar6 (TACC):

* TPV36 and TPV37 at 50 m resolution: 4.7 hours on 512 CPUs.
* TPV104, a 15 s simulation (1875 time steps): 0.4 hours on 40 CPUs.

## Historical full-run wall clock, by case

Full-term timings (mesh spacing / MPI ranks / seconds for fortran, numpy,
jax) recorded before the test suite moved to a single, shorter comparison
term. Useful as a rough guide to relative cost between cases, not as a
live benchmark:

| case | dx (m) | ranks | fortran (s) | numpy (s) | jax (s) |
|---|---|---|---|---|---|
| test.tpv8 | 500 | 4 | 13.3 | 74.3 | 20.4 |
| test.tpv10 | 500 | 4 | 39.0 | 331.3 | 98.5 |
| test.tpv104 | 500 | 4 | 33.6 | 309.5 | 86.1 |
| test.tpv1053d | 500 | 4 | 51.0 | 368.3 | 98.1 |
| test.drv.a6 | 500 | 4 | 96.5 | 687.6 | 209.5 |
| test.meng2023a | 400 | 4 | 49.7 | 400.1 | 99.0 |
| test.meng2023cb | 400 | 4 | 50.7 | 346.3 | 113.9 |
| test.tpv29 | 500 | 4 | 148.3 | 1108.5 | 229.4 |

## Multi-fault python-jax-mpi

Measured 2026-10-03 (`docs/perf_ledger.jsonl`), 4 ranks, real MPI (one
process per rank), `test.tpv22`/`test.tpv23` (two-fault, their own 15 s
term):

| case | ranks | wall (s) |
|---|---|---|
| test.tpv22 | 4 | 585-589 (contended box, busy_ceiling 0.5) |
| test.tpv23 | 4 | 276-280 (contended box, busy_ceiling 0.5) |

Both run under the same `MPI4NodalQuant.DECOMP` box decomposition as every
other python-jax-mpi case; multi-fault adds a per-row divide mask and a
per-fault master-id lookup (`src/python/eqdyna/MPI4NodalQuant.py`) rather
than a different parallelization scheme, so no separate scaling curve is
expected or measured here -- see [Core scaling](#core-scaling) above for the
single-fault curve this case set shares.

## HPC scaling suite

`testsys/perf/hpc_scaling_suite.py` builds, submits, collects and analyzes a
strong/weak/size scaling sweep on an HPC machine registered in
`scripts/machines.py`:

```
python3 testsys/perf/hpc_scaling_suite.py --submit --machine ls6 --account <acct> \
    --backend fortran --test strong --level medium
python3 testsys/perf/hpc_scaling_suite.py --collect --suite-dir <dir from --submit>
python3 testsys/perf/hpc_scaling_suite.py --analyze <tarball from --collect>
```

`--test strong` holds total element count fixed while ranks double (only
`par.nx/ny/nz` change, not `par.dx`); `--test weak` holds elements-per-core
fixed as the rank grid grows (`par.dx` refines); `--test size` holds core
count fixed and sweeps `par.dx` to find the elements-per-core/bytes-per-element
sweet spot before memory runs out. `--level big` (10^10 elements) refuses
`--backend python-jax-mpi` in code -- Fortran only at that size. `--analyze`
reports each rank's own per-step cost (via the one shared per-step-by-difference
routine, `testsys/perf/run_numa_scaling.per_step_and_fixed`), a per-step bucket
breakdown (element/fault/exchange/wait), the setup+io fixed cost read directly
off the longer run (never differenced -- a one-time cost does not scale with
step count), max/mean imbalance across ranks, the 128->256-node step boundary
when present, and a parity check of the first strong-scaling point against its
registered e2e reference via `testsys/compare.py`.

The actual until-OOM `--test size` sweep and any real HPC submission are left
to real HPC hardware; this box is memory-constrained and does not run them.
`testsys/regression/test_hpc_scaling_suite.py` unit-tests the collect/analyze
arithmetic against a committed profile fixture instead.

## Hardware assumptions

The figures above were measured on a 64-core (2 sockets x 32 cores, 8 NUMA
nodes of 8 cores each) Linux workstation. Runtime on other hardware will
differ, especially for MPI scaling and for GPU runs of the JAX backend.

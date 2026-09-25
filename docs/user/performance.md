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

Running every benchmark case on every backend (`python3 testsys/run.py all`)
on a 64-core workstation took between about 2100 and 2300 seconds across
three runs between 2026-09-17 and 2026-09-19. Wall clock depends heavily on
other work sharing the same machine at the time. Cases are independent and
each uses roughly one core, so throughput comes from running many cases at
once, not from any single case scaling across cores.

## Core scaling

Core-count scaling of a single JAX-CPU run is not currently well
characterised: one measurement on `test.tpv8` shows about 14x from one
unpinned core to several, but an earlier, more detailed scaling curve
predates a solver restructure and has not been reproduced. Treat any core
scaling claim as provisional until it is remeasured on current code.

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

## Hardware assumptions

The figures above were measured on a 64-core (2 sockets x 32 cores, 8 NUMA
nodes of 8 cores each) Linux workstation. Runtime on other hardware will
differ, especially for MPI scaling and for GPU runs of the JAX backend.

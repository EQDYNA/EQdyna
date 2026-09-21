# item 43 -- explicit parallelism for the jax backend (worktree wt-jaxshard, branch jax-shardmap)

Baseline being beaten (pathway item 43, 2026-09-19, v5.12.0/737358f, test.tpv104,
231040 cells, friclaw 4, per-step by n_lo/n_hi difference, compact numactl pin):
Fortran 931.65/466.64/237.62/117.07/64.74 ms/step at 1/2/4/8/16 -> 14.39x.
jax 611.00/410.27/335.19/253.79/221.55 -> 2.76x.

## H1 -- jax-native explicit sharding (shard_map over host CPU devices). DONE, bounded.
backend.run_time_loop_sharded + driver's nodal_sync at driver.f90:27's position.
Element arrays cut on their leading axis; nodal stages replicated; one psum/step.
- Works. Physics moves only by reassociation (tpv8, 20 steps, 2 devices vs 1:
  max abs diff 2.7e-09 in one frt column, <= 1e-15 in the rest).
- The duplicate-index scatter-add DOES shard explicitly. Item 43's "0 of 38
  scatters annotated" was about XLA's AUTOMATIC partitioner and does NOT bind
  explicit sharding. That question is now answered; it was worth asking.
- MEASURED collective cost, pure psum, no other work, 3.82e6 float64 = 30.5 MB
  (scratch microbenchmark, numactl-pinned, N cpus for N devices):
    N=4 14.8 ms | N=8 18.0 ms | N=16 19.4 ms
  Small messages are cheap: 1 MB costs 1.4 ms at 16 devices, 8 doubles 0.21 ms.
- Same-session sweep (paused when the MPI redirect arrived): element dev=1
  895.98, dev=2 424.14 ms/step (2.11x). dev=1 is 47% slower than the 09-19
  611.00 on the same metric -- today's tenancy, which is exactly why a same-day
  Fortran denominator is mandatory.
- VERDICT: bounded by two things a device mesh cannot fix -- the nodal stages
  stay replicated on every device (Amdahl), and the collective is O(NEQ) where
  MPI moves O(boundary). Kept as the comparison record, per owner ruling.

## H2 -- real MPI, one process per rank, jax owns the local element kernel. IN PROGRESS.
src/python/eqdyna/MPI4NodalQuant.py (decompose + exchange) + driver.run_mpi.
Design: every rank builds the PROVEN serial mesh and RESTRICTS it (contiguous
work-weighted element slab; nodal rows = the nodes its elements touch; the fault
rows it owns), so all index arrays come from the gated serial path and each rank
writes frt.txt<rank>, which testsys/frt_canonical.py already globs.
- The step is jitted in TWO halves either side of driver.f90:27; MPI4NodalQuant
  runs between them on the host. The halo gather is an extra OUTPUT of part_a and
  the received delta an extra ARGUMENT of part_b, so per-step host traffic is
  O(boundary) and the 30.5 MB force array never crosses to the host.
- Two sync modes, both correct, both measured: halo (Sendrecv, ascending rank,
  deterministic) and allreduce (O(NEQ), no ownership bookkeeping at all).
- HALO SIZE IS SMALL, measured not assumed: test.tpv8 at 2 ranks, 13413 of
  719643 equations = 1.86%. Serial element order is spatially coherent, so a
  contiguous slab has a surface and not a volume.
- Smoke: 1/2/4 ranks run; per-rank frt files written; rank and element counts
  print for every rank.
- TRAP FOUND AND FIXED: jax dispatch is asynchronous, so timing the exchange
  without block_until_ready charged the whole element kernel to MPI (296 of 622
  ms/step on a ONE-rank run with zero neighbours). The barrier is now timed
  separately from the exchange, so load imbalance cannot be read as collective
  cost.
- TRAP FOUND AND FIXED: breaking the fused fori_loop open costs 3x per step
  unless the carry is DONATED (11 arrays, ~120 MB of copy per step on tpv104).
- OPEN: MPI_Allreduce reduction order is not guaranteed reproducible. Must be
  checked bit-for-bit before that mode can be a gated path.

## Box conditions (they matter here more than the code does)
- Thread oversubscription is NOT the problem: pinned to one cpu a rank holds
  6-10 threads but EFFECTIVE_CORES reads 0.96-1.00. Unpinned (8 cpus/rank,
  mpirun default) it is 34 threads.
- --bind-to core --map-by core is WORSE than nothing: it always takes socket-0
  cores 0..N-1, some of which carry a tenant at 1.00 utilisation, and a
  halo-synchronised solver runs at the speed of its slowest rank. Placement is
  now --cpu-set <k least-loaded cpus, measured now> --bind-to cpu-list:ordered,
  with each chosen cpu's busy fraction printed beside the number it produced.
- A strict 0.20 busy ceiling selects nothing today (0 of 64 cpus, twice).

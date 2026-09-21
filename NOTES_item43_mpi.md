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
- ANSWERED 09-21: allreduce reproducibility. test.tpv8, 4 ranks, 20 steps, two
  runs per mode, compared TWO ways:
    * frt.txt0..3 byte-for-byte: identical a vs b in BOTH modes, and identical
      ACROSS modes (allreduce == halo).
    * float64 state (sha256 of velArr/dispArr/force/fric/fnft as raw bytes,
      because frt is E18.7E4 = 7 digits and a byte-equal frt does NOT prove the
      reduction was bit-stable): allreduce identical a vs b, 4/4 ranks hashed.
  VERDICT: bit-stable as measured. But MPI does not PROMISE reduction order,
  and this is one box, one Open MPI build, one message size. So: the gated
  path is sync=halo, which is deterministic BY CONSTRUCTION (Sendrecv with
  each neighbour in ascending rank order, summed in that order); allreduce
  stays a measurement mode. The ordering guarantee is a property of the code
  in halo mode and only an observation in allreduce mode.
  Probe kept at scratchpad hashrun.py; the frt-level half is what the gate
  would run, so it belongs in the MPI cell's own test, not here.
- MEASUREMENT TRAP (third of the class): `grep '^HASH'` dropped one rank of
  four because the X11 'Invalid MIT-MAGIC-COOKIE-1 key' noise mpirun emits has
  no trailing newline and prefixes a rank's line. Two empty files then diffed
  EQUAL -- a vacuously green reproducibility check. The probe now counts the
  hashes and fails if it did not get one per rank. Same shape as the -np
  guard already in run_mpi_scaling.py: a multi-process comparison goes vacuous
  by losing ranks, not by erroring.
- RUNNER BUG FIXED: eqdyna3d's per-rank summary line did not print
  wait_ms_per_step or sync, so run_mpi_scaling.py died with KeyError after the
  first 160-step point (log mpi_scaling_2026-09-21_tpv104.log, sha c4afa78).
  Both keys are in the print now.

## The gate gap, and the invocation contract a cell would need
The 30-cell sweep CANNOT exercise this path as built: testsys/e2e/run_e2e.py's
`run_cell` sends every python backend through `run_standalone`, which is
literally `python3 -m eqdyna <case_dir> --backend <numpy|jax>` in-process --
one rank, no mpirun, by construction. So no amount of passing the existing
sweep says anything about driver.run_mpi. What the cell needs (I am NOT
building it; this is the contract for whoever does):

1. A FOURTH value on the backend axis, `python-jax-mpi`, not a flag on
   `python-jax`. Backend is the sweep's axis and a cell is (case, backend);
   a hidden flag would make one cell name mean two different solves.
2. Invocation:
     mpirun --bind-to none -np <matrix.PY_MPI_RANKS[case]> \
       python3 -m eqdyna <case_dir> --backend jax --mpi
   with env EQDYNA_MPI_SYNC=halo. `halo` is the deterministic-by-construction
   mode (see the reproducibility entry above); the gate must PIN it rather
   than inherit whatever the environment holds -- a gated path that reads a
   mode from the environment is two paths.
3. Case build: `make_serial_case` unchanged. Its "forced serial" edit writes
   Fortran decomposition parameters, which the port does not read -- every
   rank builds the full serial mesh and restricts it. Worth one assertion in
   the cell that the python ranks' summed element count equals the serial
   totalNumOfElements, so a decomposition that quietly drops elements cannot
   pass by writing a short frt.
4. Comparison: NONE OF IT IS NEW. Each rank writes `frt.txt<rank>` holding the
   nodes it OWNS, `frt_canonical.frt_rank_files` already globs `frt.txt*` and
   sorts numerically, and `compare.compare_cell` already compares the
   canonical array against the one committed reference at that case's
   CASE_BOUND. A 4-rank python run compares like a 4-rank Fortran run.
5. The vacuity guard the cell must carry: assert it counted
   `PY_MPI_RANKS[case]` frt files AND that the ranks together owned
   `nftnd` fault nodes. driver.run_mpi already raises on the second (it
   allreduces the owned count against nftnd), so the cell needs the first.
   Without it, a launch that silently started 1 rank produces a perfectly
   green canonical comparison -- because a 1-rank run of this path is also
   correct, just not what the cell claims to test.
6. Cost: it is 10 new cells if added for every case, and CI memory is the
   binding constraint (the python cells were admitted at 2.91-4.34 GB
   EACH, times the rank count here). Start with the cheap cases and
   register in matrix.py's GATE, not CI_CELLS, until measured.

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

# Gate evidence for `01a1040` (jax-MPI rank-local subdomain)

**WHAT:** the full-suite run that gated `01a1040` before it landed on master,
plus the structured per-cell record from the same run.

**WHY IT IS HERE:** this evidence existed only in a session scratchpad, which
does not survive the session. Item 45 and the clean-room v5.2.0 baseline are on
the board because that failure mode has already nearly cost this project
irreplaceable data. A gate that cannot be re-read is a gate you have to take on
trust.

**FILES**
- `run_all_sweep_2026-09-22.log.gz` — `python3 testsys/run.py all`, gzipped
  (2,289,425 bytes raw, md5 `9fb4b439864c5836a69d4df8e3af4edf`; verified to
  round-trip byte-for-byte before committing). 31 of 40 cells ran, 31 passed,
  0 failed, 9 declared-unsupported, unit/regression/e2e all SUCCESS, 2904.6 s.
- `e2e_cells_2026-09-22_093928_2161337.json` — the per-cell wall-clock snapshot
  the harness wrote from that same run; its rows are also in
  `docs/perf_ledger.jsonl`.

**THE CELL THIS GATE EXISTS FOR:** `test.tpv8 x python-jax-mpi`,
`max|diff| = 1.220703e-10` against bound `1.0e-08` — unchanged to every printed
digit by a 3,250-line index remap. Per rule 5a that bound was not the operative
gate: the landing was held to BIT-IDENTITY, all four `frt.txt<rank>` md5-equal
against master's pre-landing solver on the same inputs and placement.

**VERIFY**
    gunzip -c run_all_sweep_2026-09-22.log.gz | md5sum   # 9fb4b439864c5836a69d4df8e3af4edf
    gunzip -c run_all_sweep_2026-09-22.log.gz | tail -20 # the SUMMARY block

# Gate evidence for `a16cc98` (item 58: one GPU cell per card)

**WHAT:** the full-suite run that gated `a16cc98`, plus the GPU before/after
transcripts — which are the part a CPU sweep cannot show.

**FILES**
- `run_all_sweep_2026-09-22.log.gz` — `python3 testsys/run.py all` on the CPU
  path (2,292,364 bytes raw, md5 `e8cfe9fa720ff7103dc72e1504807974`, round-trip
  verified). 31 of 40 cells ran, 31 passed, 0 failed, unit/regression/e2e all
  SUCCESS, 3340.2 s.
- `gpu_before_four_on_one_card.txt` — the FAILURE, reproduced before the fix:
  four concurrent raw `python3 -m eqdyna --backend jax` under
  `JAX_PLATFORMS=cuda`, `CUDA_VISIBLE_DEVICES=0` and no explicit memory
  fraction. Two of the four died, one with `INTERNAL: Failed to allocate
  11406984 bytes for new constant`, one with `RESOURCE_EXHAUSTED: Failed to
  load in-memory CUBIN ... CUDA_ERROR_OUT_OF_MEMORY`.
- `gpu_after_one_card.log.gz` — the same card, after the fix: both cells
  forced onto device 0 with `EQDYNA_E2E_GPUS=0`, serialised by the slot pool.
  `test.tpv8` SUCCESS 11.3 s at `max|diff| = 4.425048e-10` (bound 1.0e-08),
  `test.tpv36` SUCCESS 62.6 s at `5.359716e-09` (bound 1.0e-06).

**THE BOUND ON THE CLAIM.** With only TWO concurrent processes on one card the
failure does NOT reproduce — both exited 0. XLA's 0.75 is of FREE memory, so
on a 40 GB A100 the second process still fits. The collision needs more than
two cells sharing a card, which is what a 31-cell sweep produces and what the
2026-09-22 GPU sweep hit. Recorded so nobody later reads "two cells cannot
share a card" as the mechanism.

**WHAT WAS NOT TOUCHED:** `test.tpv36` itself — no case file, no reference, no
bound. The failure was never parity: run alone the cell passed at 7.06e-09
against 1.0e-06, and after this fix it passes at 5.359716e-09 on the very card
that OOM'd.

**VERIFY**
    gunzip -c run_all_sweep_2026-09-22.log.gz | md5sum   # e8cfe9fa720ff7103dc72e1504807974
    gunzip -c run_all_sweep_2026-09-22.log.gz | tail -20

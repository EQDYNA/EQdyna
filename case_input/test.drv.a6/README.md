# test.drv.a6

Internal case (no SCEC spec): a fractal-rough strike-slip fault with
Drucker-Prager off-fault viscoplasticity and PML absorbing boundaries,
used for deterministic ground-motion studies. Verifies the plasticity +
PML + rough-fault path (`pastReleaseNotes.md`: "a fractal fault plastic
model for ground motion application").


> **WITHDRAWN FROM THE SWEEP, 2026-09-15 — reference kept.** This case is
> `matrix.REFERENCE_ONLY`: it is no longer gated, and a green e2e sweep says
> nothing about it. The reason is a numerics question, not a tolerance
> question — Fortran against its own committed 4-rank reference, changing
> nothing but the decomposition, gives **329 arrival flips of 5151** and a
> **2.12e7** max difference (p50 4.03e-03; median arrival shift 0.0417 s =
> one time step). A case whose Fortran-vs-Fortran floor is that high cannot
> certify a backend. Everything stays in place for its return: this
> reference, `DRV_A6`, `compare.flip_budget_gate`, and
> `testsys/parity/evidence_drv_a6_chaos.py`. The two experiments that decide
> it — a 1/2/4/8-rank comparison, and the decisive dx=250 refinement — are
> written up as item 32 in `pathway_forward.md`.

![reference result at fast-tier resolution](cRuptureDynamics.png)

| tier | dx (m) | term (s) | ranks (decomp) | ~cells | wall time | reference |
|---|---|---|---|---|---|---|
| fast (was the gate) | 500 | 5 | 4 (2,2,1) | 256k | ~48s (suite-aggregate, not case-isolated) | reference kept at `test.reference.results/test.drv.a6/`, NOT gated since 2026-09-15 — the chaos-aware gate it used (`testsys/matrix.py`'s `DRV_A6`, `testsys/compare.py`'s `flip_budget_gate`) is still implemented and still regenerable by `evidence_drv_a6_chaos.py`, just not wired into the sweep |
| full | -- | -- | -- | -- | -- | EXCLUDED: no published full-resolution run in README.md/pastReleaseNotes.md beyond this test config |

Provenance: SHA 5e76e8e, cotopaxi, Ubuntu 22.04/gfortran 11.4.0/OpenMPI
4.1.1, 2026-09-14, shared box (fast wall time is a 7-case suite total /7,
not individually timed).

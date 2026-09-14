# test.meng2023a

Internal case (no SCEC spec): a strike-slip fault in a layered 1D velocity
structure (a Cushing-earthquake-style material profile), slip-weakening
friction, single central nucleation patch. Verifies heterogeneous-material
(`par.mat`) handling (`pastReleaseNotes.md`: "test system now supports ...
meng2023a, meng2023cb").

| tier | dx (m) | term (s) | ranks (decomp) | ~cells | wall time | reference |
|---|---|---|---|---|---|---|
| fast (gate) | 400 | 5 | 4 (2,2,1) | 510k | ~48s (suite-aggregate, not case-isolated) | frozen golden, `test.reference.results/test.meng2023a/` |
| full | -- | -- | -- | -- | -- | EXCLUDED: no published full-resolution run in README.md/pastReleaseNotes.md beyond this test config |

Provenance: SHA 5e76e8e, cotopaxi, Ubuntu 22.04/gfortran 11.4.0/OpenMPI
4.1.1, 2026-09-14, shared box (fast wall time is a 7-case suite total /7,
not individually timed).

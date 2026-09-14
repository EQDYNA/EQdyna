# test.tpv1053d

SCEC TPV105-3D: a vertical strike-slip fault with rate-and-state strong
velocity-weakening friction plus thermal pressurization (shear heating
raises pore pressure, weakening the fault), verifies the TP module against
the cross-code SCEC reference solutions.

| tier | dx (m) | term (s) | ranks (decomp) | ~cells | wall time | reference |
|---|---|---|---|---|---|---|
| fast (gate) | 500 | 5 | 4 (2,2,1) | 282k | ~48s (suite-aggregate, not case-isolated) | frozen golden, `test.reference.results/test.tpv1053d/` |
| full (spec) | -- | -- | -- | -- | -- | EXCLUDED (see below) |

Full tier excluded: the SCEC description+file-formats PDFs
(https://strike.scec.org/cvws/tpv105_3D_docs.html) state the 0-15 s
duration but never a required element size -- only a free-text field for
submitters to report their own. No number to run at without inventing one.
Cross-code results: https://strike.scec.org/cvws/metric_cvv1_u1/tpv105-3d/metric_cvv1_tpv105-3d_ar_0.html

Provenance: SHA 5e76e8e, cotopaxi, Ubuntu 22.04/gfortran 11.4.0/OpenMPI
4.1.1, 2026-09-14, shared box (fast wall time is a 7-case suite total /7,
not individually timed).

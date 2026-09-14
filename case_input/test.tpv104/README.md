# test.tpv104

SCEC TPV104: a vertical strike-slip fault in a homogeneous elastic
halfspace with rate-and-state (ageing law) friction, verifies the RSF
implementation (state-variable evolution, velocity-dependent friction)
against the cross-code SCEC reference solutions.

| tier | dx (m) | term (s) | ranks (decomp) | ~cells | wall time | reference |
|---|---|---|---|---|---|---|
| fast (gate) | 500 | 5 | 4 (2,2,1) | 231k | ~48s (suite-aggregate, not case-isolated) | frozen golden, `test.reference.results/test.tpv104/` |
| full (spec) | 50 | 12 | 16 (4,2,2) | ~231M (est.) | ~80-120h (est.) | report-only, no golden reference at this resolution |

Spec: recommended element-size/node-spacing 50 m, time histories 0-12 s
(https://strike.scec.org/cvws/tpv103_104docs.html, download/SCEC_validation_slip_law.pdf).
Cross-code results: https://strike.scec.org/cvws/metric_cvv1_u1/tpv104/metric_cvv1_tpv104_ar_0.html

Provenance: SHA 5e76e8e, cotopaxi, Ubuntu 22.04/gfortran 11.4.0/OpenMPI
4.1.1, 2026-09-14, shared box (fast wall time is a 7-case suite total /7,
not individually timed); full-tier wall time is `(500/dx)^4 * term-ratio`
scaled from that aggregate -- an estimate, not a measurement, and the
50 m/12 s spec is ~13 days extrapolated, worth an owner sanity-check
before actually launching.

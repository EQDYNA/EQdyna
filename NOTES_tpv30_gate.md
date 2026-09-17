# TPV30 promotion — working notes (rule 17)

Worktree: agent-a0f87ff963afe7929. Base: master 03ba055 (matches origin).

## Scoping findings (read before trusting anything below)

- Item 24(a) (G6 "half-traction BLOCKER" in the old scratch draft) is
  ALREADY RESOLVED as of commit c7c4f5f (2026-09-16): measured on drv.a6,
  Tn ratio 1.0018, Ts ratio 0.9858 -- NOT half. Root cause was item 26's
  `arn` double-count, fixed 2026-09-14. The old draft's G6 note is STALE;
  it predates that fix. Not re-litigated here; just confirmed via git log
  and will re-verify empirically for TPV30 itself once a run exists.
- Item 24(b)(c)(e)(f) landed fa180b4 (2026-09-17): par.viscoplasticRelaxTime,
  par.devStrTaperDepthStart/End, lib.shearModulusFromPar, par.plasticOutputHalfWidth
  are all real case inputs now. G1/G2 (old draft) are closed at the src level.
- G3 (+7.3215 m constant): untouched by owner decision (item 24(d)), LOAD-BEARING
  per c7c4f5f's own measurement -- do not touch.
- G5 (b11+b33 = 1.999999 not 2): negligible (164 Pa @ 10km), not addressed by
  the spec itself, carried forward as-is.
- G4 (plastic-strain output window box) now a parameter (par.plasticOutputHalfWidth)
  -- TPV30 compset must set it explicitly (40x20 km fault >> the 5x2x8 km default).
- Spec confirmed (TPV29_30_Description_v06, read directly, pages 11/14-18):
  p.11 "The material properties are the only difference between the elastic
  benchmark (TPV29), and the viscoplastic benchmark (TPV30)." p.16 "Benchmarks
  TPV29 and TPV30 use linear slip-weakening friction" (same 6 params, same T(r)
  nucleation formula, Part 5/6 apply to BOTH benchmarks identically) -- so
  TPV==30 belongs in swtwNucleation's branch list exactly as TPV==29 does.
  p.11: cohesion c=1.18 MPa, bulk friction nu=0.1680, Tv=0.05s. p.13: Omega(depth)
  taper 17000-22000 m. All match the old draft's own transcription.
- faulting.f90:405 and faulting.py:109 confirmed missing TPV==30 (grepped).
- tpv29GeometryTools.py (case_input/test.tpv29/) is fully generic -- no
  TPV29-specific hardcoding beyond file names/docstrings. Reusable verbatim
  for TPV30 (same official rough surface, spec-confirmed).
- case.setup / scripts/lib.py's ensureFaultRoughGeometryForCase /
  requireFaultGeometryResolution are generic, no changes needed there.

## Plan (rule 17 order)

1. [x] Spec fetch -- already satisfied (scratch/tpv29/downloads/TPV29_30_Description_v06.pdf),
   confirmed by direct read, cited above.
2. [ ] Geometry: case_input/test.tpv30/ ships tpv29GeometryTools.py (byte copy)
   + bFault_Rough_Geometry.tpv29.100m.txt (byte copy) -- same official surface.
   50m file NOT shipped yet (13.9 MB, not needed until the full/spec tier is
   actually run -- rule 17 step 5 says record, not run).
3. [ ] Code branch: add `.or. TPV==30` to faulting.f90:405 AND the tuple in
   faulting.py:109. par.tpv = 30 in the real case (not 36).
4. [ ] Gate: dx=500, 4 ranks (2,1,2 like tpv29), term=20. Freeze frt.canonical.txt.
   Register testNameList.py + testsys/matrix.py (CASE_BOUND measured after a
   real run, GATE='abs-max').
5. [ ] full_specs.py: dx=50 entry, citation to Part 3 p.15 "50 m node spacing...
   100 m node spacing" -- record only, do not run.
6. [ ] Promote scratch/tpv30/compareTpv29Tpv30.py -> testsys/parity/
   evidence_tpv30_vs_tpv29_contrast.py, report-only pattern (like
   evidence_tpv29_scec_comparison.py: never asserts, exits 0, prints
   REPRODUCES/DRIFTED-style verdict, writes gitignored JSON).
7. [ ] Full sweep, all 3 backends, test.tpv30 must be SUPPORTED on every one
   (no UNSUPPORTED entries for a new case, per rule 17's own text).

STATUS UPDATE (2026-09-17, mid-session):

- [x] Step 2: case_input/test.tpv30/ created -- tpv29GeometryTools.py and
  bFault_Rough_Geometry.tpv29.100m.txt byte-copied (md5
  f6df5ececc9d95d0c6db43c9ba0c3c6e, matches test.tpv29's). user_defined_params.py
  written: par.dx=500 fast gate, par.tpv=30, viscoplasticRelaxTime=0.05,
  devStrTaperDepthStart/End=17000/22000, plasticOutputHalfWidth=(25e3,20e3,22e3)
  (fixes G4 for this fault size). Lazy faultGeometryWriter pattern (NOT
  write-at-import like the old scratch draft). README.md written.
- [x] Step 3: `.or. TPV==30` added to faulting.f90:405 AND the tuple in
  faulting.py:109 (was `(29,36,37,201)`, now `(29,30,36,37,201)`).
- [x] Step 4 (partial): built fresh (`./install-eqdyna.sh -m ubuntu`), ran
  test.tpv30 manually first (create.newcase, case.setup, mpirun -np 4 eqdyna,
  plotRuptureDynamics) to sanity check BEFORE trusting the harness:
    - case.setup printed VISCOPLASTIC: Tv=0.05s, taper 17000->22000m,
      plasticOutputHalfWidth as configured -- wiring confirmed live.
    - geometry validated clean (nnx=81,nnz=41,dx=500, same corner/roughness
      numbers as test.tpv29's own dx=500 gate).
    - run completed exit 0, term=20s reached.
    - frt.canonical.txt: 3321 rows, 2984/3321 nodes nucleated (TPV29's own
      reference nucleates 2974/3321 -- close, expected small difference from
      added plasticity, NOT the 0/3321 degenerate-branch bug rule 17 step 3
      warns about). Max final slip 2.91 m vs TPV29's reported ~2.90 m at
      100m/~2.78 at 500m -- physically sane, not a G6-style half-traction
      symptom.
  Froze test.reference.results/test.tpv30/frt.canonical.txt from this run.
  Registered testNameList.py (+test.tpv30, 4 ranks) and testsys/matrix.py
  (CASE_BOUND provisional 1e-6, GATE='abs-max', same family as tpv29/36/37) --
  needed so matrix.py imports at all with test.tpv30 in CASES.
- [x] Step 5: full_specs.py FULL_SPECS['test.tpv30'] added, dx=50 term=20
  nx=4 ny=1 nz=4, citation to spec Part 3 p.15. NOT run (only the 100m
  geometry file is shipped in case_input/test.tpv30/ today; noted in README).
- [x] Step 6: testsys/parity/evidence_tpv30_vs_tpv29_contrast.py written
  (report-only, never asserts, exits 0), promoted from
  scratch/tpv30/compareTpv29Tpv30.py, rewritten to read dx/fault-extent/
  hypocenter from each case's OWN compset rather than hardcoding, and to
  drop the stale G6 blocker language. NOT YET RUN (needs completed
  test/test.tpv29 + test/test.tpv30 fortran-cell dirs from the sweep below).
- [~] Step 7: `python3 testsys/e2e/run_e2e.py --cases test.tpv30
  --backends fortran,python-numpy,python-jax` LAUNCHED, running in background
  (all 3 backends concurrently). Will report actual observed diffs and
  tighten CASE_BOUND from the provisional 1e-6 once it completes, then
  re-run to confirm the tightened bound still passes (matches how
  tpv36/tpv37 bounds were set). THEN run `python3 testsys/run.py all` fresh
  for the full 33-cell table.

NEXT: read sweep output, tighten bound, re-verify, run evidence script,
run full gate, update this file again before finishing.

## STOP — finding, not a landing (2026-09-17, final)

The full sweep (`run_e2e.py --cases test.tpv30 --backends
fortran,python-numpy,python-jax`) ran to completion and ALL THREE cells
FAILED against the provisional 1e-6 abs-max bound:

    test.tpv30  fortran       FAIL (missing reference nc -- fixed by copying
                              fault.dyna.r.nc from this same run; not a real defect)
    test.tpv30  python-numpy  FAIL  max|diff|=4.035089e+08  bound=1.0e-06
    test.tpv30  python-jax    FAIL  max|diff|=4.035089e+08  bound=1.0e-06  (IDENTICAL numbers)

python-numpy and python-jax produced BIT-IDENTICAL worst-case numbers
(same row, same column, same ref/run values to all printed digits) and
their FULL canonical grids agree to 1e-3 max over 3321x22 values. Both
disagree with the Fortran reference by up to 4.0e8 Pa (~30% relative) at
3310 of 3321 nodes.

Binary search performed (rule: never conclude "just inaccurate" without
looking):
  - Built two throwaway diagnostic case dirs (test/diag_short_f,
    test/diag_short_p) at par.term=1.0 and 6.0 (everything else identical
    to the real compset) to localize WHEN the divergence starts, since the
    full 20 s run is too coarse a unit to bisect by re-running.
  - t=1.0 s (24 steps): fortran (4-rank) vs python-numpy (serial) canonical
    grids agree to 1.6e-14 -- BIT EXACT, all 3321 nodes, all 22 columns,
    including the hypocenter. This directly confirms item 24(b)/(c)'s new
    plumbing (Tv, taper, b11/b33/theta initial stress, the Drucker-Prager
    return map itself) is bit-correctly ported -- ELIMINATED as the cause.
  - t=6.0 s (144 steps): already diverged. Max diff 1.6e8 Pa; column 3
    (fnft) shows a >9.9e4 s difference at >=1 node -- a genuine
    rupture-arrival flip (same signature as test.drv.a6's known
    bistability). ALL 3321 nodes already differ by >1e3 Pa.
  - Did NOT narrow further inside the 24-144 step window (would need
    per-step dumps / a debug Fortran build on both sides -- correctly the
    next step per my own protocol, but a real time investment beyond this
    pass's budget).

CONCLUSION: this is a real, deterministic (NOT per-backend-chaotic --
numpy==jax to 1e-3 while both are 30% off Fortran rules out independent
floating-point chaos as the primary mechanism) divergence in how the
ported Drucker-Prager viscoplastic return map interacts with the ROUGH
(non-planar) fault over many steps. test.drv.a6 is the only other gated
case combining C_elastic=0 with a rough fault, and it needed a dedicated
flip-budget gate (DRV_A6 dict in matrix.py, built from real measurement)
rather than abs-max -- TPV30 may need the same, or a real fix. Neither
attempted here; per the mission brief's own stop-trigger ("if the scope is
bigger than this brief assumes... STOP and report"), this qualifies:
the brief assumed TPV30 could reuse TPV29/36/37's abs-max gate pattern
once G6 was confirmed resolved (it is) and Tv/taper were wired (they are,
verified bit-exact) -- it did not anticipate a fresh algorithmic
divergence surfacing only when rough-fault geometry and viscoplasticity
are combined for the first time under friclaw=1.

ACTION TAKEN: reverted the sweep registration only.
  - testNameList.py: test.tpv30 REMOVED (with a comment explaining why and
    pointing here).
  - testsys/matrix.py: CASE_BOUND/GATE entries for test.tpv30 REMOVED (with
    the same explanation, plus the concrete measured numbers).
  - `python3 testsys/run.py unit regression` re-run after the revert:
    SUCCESS both tiers -- the existing 30-cell table is untouched by this
    work (I did not have session budget to re-run the full ~1500s e2e
    sweep after the revert; the revert itself is two clean line-removals in
    files nobody else's cells depend on, reviewable directly).

KEPT (real, verified progress, not reverted):
  - src/fortran/faulting.f90 / src/python/eqdyna/faulting.py: `TPV==30`
    added to swtwNucleation's branch list on BOTH sides (rule 17 step 3).
    Verified independently correct via the spec PDF (p.16: "Benchmarks
    TPV29 and TPV30 use linear slip-weakening friction") and via the
    t=1s bit-exact diagnostic above (which exercises this formula).
  - case_input/test.tpv30/ (compset, geometry, README with the finding
    above documented as the compset's own gate-status section).
  - test.reference.results/test.tpv30/ (frt.canonical.txt, fault.dyna.r.nc,
    cRuptureDynamics.png) -- a verified-correct Fortran run, kept as a
    candidate reference for whoever closes the finding, NOT wired into the
    gate.
  - testsys/e2e/full_specs.py FULL_SPECS['test.tpv30'] (rule 17 step 5,
    record-only, independent of gate status).
  - testsys/parity/evidence_tpv30_vs_tpv29_contrast.py (rule 17 step 6,
    report-only; not run to completion against a real pair of dirs as part
    of this finding, but syntactically verified and ready).

NOT reused/left over: test/diag_short_f, test/diag_short_p, test/test.tpv30,
test/test.tpv30.python-numpy, test/test.tpv30.python-jax under the
gitignored test/ tree -- left on disk as evidence (rule 8), not committed
(test/ is gitignored).

HANDOFF: this needs either (a) a Fortran-side debug build with per-step
stress dumps to find the exact divergence, matched by an equivalent Python
dump (my own protocol's "build a debug version of the reference" step,
not done here), or (b) a deliberate flip-budget gate design measured the
way test.drv.a6's was (multiple decompositions, serial-vs-parallel Fortran
baseline, etc.) -- either is a real, separate piece of work, not a
guess-and-loosen-the-bound fix.

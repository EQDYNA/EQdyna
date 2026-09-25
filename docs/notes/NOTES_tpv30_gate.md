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

## Rule 17 step 6 satisfied: independent validation against the owner's OWN
## published TPV30 submissions (2026-09-17, follow-up mission)

New information changed what was possible: `scec_archive/tpv30/eqdyna-v3.1-
{100,50,25}m-2015/` (the owner's own 2015 SCEC/USGS cvws submissions,
gitignored, present on this box) were not available/considered when the STOP
above was written. Rule 17 step 6 ("validate against something independent,
as a committed script") is satisfiable today with data already on disk, so it
was done, WITHOUT touching the unresolved Fortran-vs-Python finding above.

**Script**: `testsys/parity/evidence_tpv30_scec_comparison.py` (new, report-
only, never asserts, exits 0 -- same pattern as
`evidence_tpv29_scec_comparison.py` and `evidence_tpv30_vs_tpv29_contrast.py`
in the same directory, which it imports from rather than duplicating generic
archive I/O, per rule 1). Run:
`python3 testsys/parity/evidence_tpv30_scec_comparison.py`.

**What it compares** (see the script's own docstring for full detail): the
ALREADY-COMMITTED `test.reference.results/test.tpv30/frt.canonical.txt`
(EQdyna Fortran, dx=500 m, the rule-17-step-4 gate run frozen earlier in this
document) against `scec_archive/tpv30/eqdyna-v3.1-100m-2015/` (the owner's OWN
2015 EQdyna submission, dx=100 m). This is Fortran-vs-Fortran, 21 years apart
-- it says NOTHING about the numpy/jax divergence above; that finding is
untouched and still open.

**Framing, stated explicitly per the mission's own bar** ("roughly match",
never an accuracy claim): the owner's own 100 m TPV30 submission itself ranks
14/14 (farthest from the group median) among 14 independent cross-code
submissions at the SCEC portal (scec_archive/tpv30/eqdyna-v3.1-100m-2015/
PROVENANCE.md), improving monotonically 100m(14/14) -> 50m(12/14) -> 25m(8/14)
in the owner's own record. So a 500 m EQdyna run matching the owner's OWN
100 m run is a REGRESSION check ("recognizably the same physics"), not an
accuracy claim, exactly the tpv36/tpv37 500 m gate framing.

**Numbers actually computed by the script (not from any README/prior claim)**:

  - METRIC 1 (rupture time, FULL 3321-node 500 m grid, exact-coordinate
    lookup into the 100 m archive grid -- dx=500 divides dx=100 exactly, no
    interpolation): 2940/3321 nodes ruptured in BOTH runs, 314 never ruptured
    in either, 44 ruptured only at 500 m, 23 ruptured only in the archive --
    99.2% of the archive's ruptured nodes also ruptured at 500 m. Median
    |dt| = 1.077 s, mean 1.476 s, max 9.611 s (coarse-grid timing noise, not
    zero, as expected).
  - METRIC 2 (ruptured area, each grid at its OWN resolution, corner rule +
    naive rule, both reused verbatim from evidence_tpv29_scec_comparison.py):
    corner rule 713.8 km2 (500 m) vs 729.9 km2 (100 m archive), ratio 0.978;
    naive rule 746.0 km2 vs 736.1 km2, ratio 1.013. (A bug was caught and
    fixed while writing this: `ref['dip'] = -ref['z']` runs DESCENDING
    (20000->0) because `ref['z']` itself runs ascending -20000->0, and
    `rupture_area_km2`'s `dz = np.diff(z).mean()` is signed -- this returned
    a NEGATIVE area on the first real run of the script, caught immediately
    because the number was obviously wrong, not silently accepted. Fixed by
    flipping to ascending order before that one call; METRIC 1/3 are
    direction-agnostic and needed no fix.)
  - METRIC 3 (final slip, all 24 archive on-fault stations, vs the 500 m
    grid's NEAREST node -- most named stations are not on 500 m multiples,
    max offset 283 m, every offset reported alongside its number): median
    percentage diff 6.3%. 20 of 24 stations are under 20%; two outliers
    (faultst167dp105 at 37.7%, faultst170dp045 at 99.3%) are both small
    absolute-slip, near-edge stations where percentage error is naturally
    amplified by a small denominator at coarse resolution -- not treated as
    evidence against the match, stated as what it is.
  - METRIC 4 (Mw, current 500 m run only -- the archive has no full 2D slip
    field, same documented limitation as TPV29's own script): M0 =
    3.876e19 N*m over all 3321 nodes, Mw = 7.026, computed fresh from this
    run's own slip field (RHO=2670, VS=3464, the case's own spec-p.11
    material), not compared cross-code (rule 1 -- no guessed baseline).

**Verdict, the two SEPARATE questions kept separate (mission Step 3, and the
script's own print_report says both, side by side, every run)**:

  (A) PHYSICS VALIDITY (this script, Fortran vs Fortran, rule 17 step 6):
  **ROUGHLY MATCHES** -- 99.2% rupture-extent overlap, ruptured-area ratio
  0.98 (corner rule), median station slip diff 6.3%. This is real and worth
  recording on its own: EQdyna's viscoplastic TPV30 physics at a coarse 500 m
  gate resolution is recognizably the same physics as the owner's own finer
  2015 submission, not a different or broken rupture.

  (B) PORT CORRECTNESS (numpy/jax vs Fortran, rule 17 step 7): **UNCHANGED,
  STILL NOT RESOLVED.** The "STOP" finding above this section stands exactly
  as written -- this script never touches it (it never runs the Python
  backends, only reads the already-frozen Fortran reference). (A) matching
  does NOT imply (B) is fixed; the numpy/jax vs Fortran divergence (up to 30%
  at t=20s, root cause not found, needs a debug-build per-step dump per the
  HANDOFF above) is completely independent of whether the Fortran physics
  itself is validated.

**Gate action: NONE.** Per the mission's own instruction, gating requires
BOTH (A) and (B); (A) now holds but (B) does not, so test.tpv30 stays
un-registered in `testNameList.py`/`testsys/matrix.py`, exactly as the
earlier STOP section left it. This is a complete, valuable combined finding
(validated physics, unresolved port divergence), not a substitute for one.

**Step 4 (costing, not running, a 100 m EQdyna validation run)**: NOT run.
Rough cost estimate from this case's own scaling, not guessed: the 500 m gate
(4 ranks, 81x41x?? mesh, term=20s) completed in the same
"minutes at 4 ranks" range rule 17 step 4 targets generally (not separately
re-timed here since it was not re-run this session). Going 500m -> 100m is a
5x reduction in dx, i.e. 5x in EACH of 2 horizontal fault-plane dimensions
that drive element count, and a similar factor in the off-fault volume mesh
plus a 5x smaller stable timestep (CFL) -- a full 3D explicit dynamic-rupture
run scales roughly as (element count) x (timestep count), so a naive
dx-only estimate is O(5^4) = 625x in raw work for the volume mesh (3
spatial dims + 1 time dimension all tightening by 5x), before accounting for
any rank-count increase to keep wall-clock down. TPV29's own 100 m run (cited
in evidence_tpv29_scec_comparison.py's module docstring) used 48 ranks and
took "tens of minutes"; TPV29's 50 m run (a further 2x tightening from its
own 100 m) was independently measured at 2-3.5 h on 16 ranks per that
script's docstring. Applying the same 500m->100m jump to TPV30 at a
similarly scaled-up rank count (48+ ranks) would plausibly land in the same
"tens of minutes to low hours" band TPV29's own 100 m run sits in, NOT the
2-3.5 h TPV29 50 m band (100 m, not 50 m, is what scec_archive ships a direct
comparison point for, per the mission's own reasoning) -- but this is an
extrapolation from a different case's own measured numbers, not a fresh
timing of TPV30 itself, and is explicitly NOT queued or started here, per the
mission's own Step 4 instruction ("cost out, do not necessarily run").

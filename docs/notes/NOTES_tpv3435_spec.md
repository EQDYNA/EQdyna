# TPV34 / TPV35 — rule 17 step 1 (spec fetch only)

Date: 2026-09-17. Model: sonnet (rate-limit override — see mission dispatch
note). Scope of this note: what the official SCEC spec PDFs actually say.
Step 3 (code-branch check) was pre-confirmed by the dispatcher before this
session started (`grep -n "TPV==34\|TPV==35" src/fortran/*.f90` → zero hits);
not re-verified here, out of this mission's scope. Steps 2, 4, 6, 7 are not
done and not started.

## TPV34 — Imperial Fault, Model 1

Confirmed against `scratch/specs/TPV34_desc.txt` (fetched from
`https://strike.scec.org/cvws/tpv34docs.html` →
`download/TPV34_Description_v10.pdf`), not just the summary table on
`benchmark_descriptions.html`:

- Single planar **vertical, right-lateral strike-slip fault**, 30 km long,
  15 km deep (Part 1, "The Imperial Fault").
- **3D velocity structure** taken from the real SCEC Community Velocity
  Model CVM-H around the actual Imperial Fault (Southern California),
  modified so the minimum Vs is clamped to 1400 m/s (Part 2).
- Linear elastic material, linear slip-weakening friction.
- Initial shear and normal stress are **proportional to the shear modulus**
  (depth/velocity-structure dependent, not a flat value like TPV8/TPV29).
- Nucleation: center of fault, 7.5 km depth, via an added shear-stress bump
  in a circular patch around the hypocenter (Part 3).
- Run duration: 0.0–20.0 s after nucleation (Part 3, "Running Time, Node
  Spacing, and Results").
- **Node spacing: NOT a single recommended number.** Spec text: "For TPV34,
  please select node spacing on the fault plane in the range of 25 m to
  50 m, and submit results for your selected node spacing." This is a
  submitter's-choice range, unlike TPV8/TPV29/TPV35's single "recommended"
  value — recorded as `EXCLUDED` in `full_specs.py`, not invented.

This is a heavier lift than a flat-stress TPV: it needs an actual CVM-H
velocity query (or a pre-extracted grid) as a material-property input, not
just uniform elastic constants. Whoever scopes the implementation should
expect a velocity-model ingestion step with no counterpart in the currently
supported TPV cases.

## TPV35 — Parkfield 2004 M6 Earthquake

Confirmed against `scratch/specs/TPV35_desc.txt` (fetched from
`https://strike.scec.org/cvws/tpv35docs.html` →
`download/TPV35_Description_v05.pdf`):

- Vertical, right-lateral strike-slip fault, 40 km long, 15.5 km deep
  (Part 1).
- **3D velocity structure = a 1D velocity profile on each side of the
  fault** (not a full 3D CVM query like TPV34) — minimum Vs 1100 m/s.
- Linear elastic, linear slip-weakening friction.
- Initial shear stress and yield stress **vary with position on the fault**,
  tuned so the resulting rupture resembles the real Parkfield 2004 Mw 6.0
  event (Ma, Custódio, Archuleta & Liu, 2008, JGR 113, B02301). This is a
  **validation benchmark**: results are meant to be compared not just
  across codes but against real seismic recordings.
- Nucleation: 10 km from one end of the fault, 8.1 km depth, via a
  **lowered yield stress** patch (not an added shear-stress bump like
  TPV34/TPV29) around the hypocenter — the nucleation stress is described
  as "built-in to the input data file," implying a supplied stress-field
  file rather than an analytic formula (Part 3).
- Run duration: 0.0–18.0 s after nucleation.
- Node spacing: **recommended 100 m**, optionally also 50 m (Part 3,
  "The recommended resolution for TPV35 is 100 meters. You may optionally
  also submit results for a resolution of 50 meters.") — recorded in
  `full_specs.py` FULL_SPECS with dx=100, term=18.
- Data files (not fetched, out of step-1 scope): `tpv35_data_files.zip`
  (input data, station locations, 1D velocity model) and two optional
  real-seismic-recording archives (NGA West 2, Ma et al. 2008) linked from
  `tpv35docs.html`. The data-files zip is almost certainly needed for step 2
  (it is the supplied stress-field/velocity input, not an SCEC-generated
  spec description) — flagging it here rather than fetching it, since this
  mission is spec-description-only.

## Original one-liner (mission dispatch) — confirmed correct

"TPV34 is understood to involve the Imperial Fault region under CVM-H 3D
velocity structure and TPV35 is understood to be a Parkfield 2004 M6
validation benchmark" — both confirmed against the fetched spec text above,
not just trusted.

## Rule 17 step-3-adjacent scoping: what code already exists to build on
(2026-09-17, wei-lin, read-only — `nmat`/`meshgen.f90`/`faulting.f90` grep and
read, no runs, done during a measurement window where nothing heavier was
permitted)

**Material model.** `nmat==1/n2mat==3` (homogeneous) and `nmat>1/n2mat==4`
(1D layered) both exist, mirrored in Fortran (`meshgen.f90:174-201`) and
Python (`meshgen.py:605-633,756-768`). The layered form keys ONLY on
`abs(elementCenterCoor(3))` (depth) — no side-of-fault or lateral
distinction. TPV35's spec ("1D velocity profile on each side of the fault")
does NOT fit this as-is: it needs a fault-relative pick (which side), not a
second depth axis. Bounded change (extend the lookup key), not a new
subsystem. **TPV34 needs a full CVM-H 3D velocity query/grid — grep for
`CVM`/`velocity_model`/`vmod` across `src/`, `scripts/` found nothing.** No
counterpart exists at all; this is the single biggest net-new piece either
case needs. Given this project's own precedent for TPV29 (ship a decimated,
provenance-cited static grid rather than a live external dependency, rule 17
step 2), the likely right shape is a pre-extracted, decimated CVM-H grid
file, not a live query library.

**Nucleation.** `swtwNucleation` (`faulting.f90:386-421`) already implements
a shared smoothed forced-rupture-time formula (radius + elapsed-time driven
friction-coefficient reduction) for TPV29/36/37/201/202. TPV34's "added
shear-stress bump in circular patch around hypocenter" reads as structurally
compatible with this family — plausibly a new `TPV==34` branch here, small.
**TPV35 is a different shape: "lowered yield stress... built-in to the input
data file"** implies an externally-supplied per-node value, not a formula.
Grepped for any existing from-scratch external per-node stress/nucleation
file input: none — the only external-file-driven per-node state in this
codebase is the mode=2 cycle restart netcdf path (`netcdf_read_on_fault_eqdyna_restart`),
a different code path with different semantics (continuing a cycle, not
seeding a single dynamic run). This is genuinely new, not a reuse.

**Validation against real recordings.** Grepped `testsys/`, `scripts/` for
any existing real-seismogram/NGA comparison capability: none (TPV29's
`evidence_tpv29_scec_comparison.py` compares against another CODE's output,
not recorded ground motion). Per this project's own pattern (rule 17 step 6
is a script, landed after gating, never a gate itself — see TPV29), this can
defer to a post-gate evidence script and is NOT a blocker for landing TPV35
on cross-code criteria alone.

**Net sequencing read, for whoever scopes the actual dispatch:** TPV35 is
closer on the material side (layering machinery ~80% there) but needs a
wholly new external per-node stress-field ingestion path with no precedent.
TPV34 is closer on the nucleation side (likely reuses `swtwNucleation`) but
needs the single largest net-new subsystem of the two, a 3D velocity grid.
Building the material-lookup generalization (fault-relative AND
external-grid-capable) as ONE shared piece rather than two bespoke paths is
the real architecture call here — and it touches code every gated case
depends on (`meshgen.f90`'s material dispatch), so it needs its own
correctness gate on the 8 already-gated cases before either TPV34 or TPV35
can build on it (rule 3). Not sized further here; still not enough to
dispatch blind, but the next scoping pass has concrete surfaces to read
instead of a blank page. `tpv35_data_files.zip` (velocity model + stress
field + stations) is still unfetched and is the next real step for TPV35;
a CVM-H extraction/decimation approach is the next real step for TPV34.

## Risk flagged for the next scoper (not fixed here, out of mission scope)

`testsys/e2e/run_e2e_full.py` iterates `FULL_SPECS.items()` directly and
calls `create.newcase` for every key. It is gated behind
`EQDYNA_FULL_LAUNCH=yes-hours` (never automatic, user-scheduled, not part of
`run.py all`/CI), so this is not a live gate risk, but the `test.tpv35`
entry added to `FULL_SPECS` in this session will make that command emit a
`create.newcase failed` failure line and a non-zero exit if anyone runs the
full tier before `case_input/test.tpv35` exists. This is the documented
trade-off of doing rule-17-step-5 "paperwork" ahead of steps 2–4; noting it
so it isn't a surprise.

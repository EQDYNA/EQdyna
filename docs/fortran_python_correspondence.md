# Fortran ↔ Python correspondence, all 26 `src/fortran/*.f90`

Audited 2026-09-23 by `sophia-okafor` (read-only, no code changed), against
`origin/master` `928c4d4`. Conductor's re-verification of the load-bearing rows
is in the last section, and **one headline finding of the audit is refuted
there** — read it before citing this table.

This is the record PROJECT_RULES **rule 23 clause 4** requires: every
`src/fortran/*.f90` is exactly one of **(a) counterpart**, **(b) documented
departure**, or **(c) declared absence**. No fourth state, and the forbidden
thing is the silent gap, not the departure. A folding is a legitimate END
STATE: the owner's ruling is that structure may depart from Fortran's module
boundaries where it buys measured performance or a backend constraint forces
it. This table exists to make the mapping VISIBLE, not to drive the port toward
one-to-one symmetry.

| Fortran file | Python counterpart / folded into (module:lines) / ABSENT | if ABSENT: what the port does instead | if ABSENT: what would catch it |
|---|---|---|---|
| `assembleGlobalKU.f90` | `assembleGlobalKU.py:261-333` | — | — |
| `assembleGlobalMass.f90` | `assembleGlobalMass.py` (whole); its embedded `MPI4NodalQuant` subroutine (`:58-245`) → `MPI4NodalQuant.py` | — | — |
| `calcB.f90` | FOLDED — `assembleGlobalKU.py:287-292` (interior), `:372-470` (`_pml`). B is never materialised; `_c()` contracts ∂N against velocity directly in Voigt order | — | — |
| `calcElemKU.f90` | FOLDED — `assembleGlobalKU.py:261-333` + `_pml:372-470` | — | — |
| `calcElemMass.f90` | FOLDED — `assembleGlobalKU.py:325-327`, gravity body force only; the `al = rdampm*velArr` inertial term is not carried | — | — |
| `calcGlobalShapeFunc.f90` | FOLDED — `assembleGlobalMass.py:126-144`, `:147-172` (`compute_element_det`), `:405-547` (`compute_element_shape`) | — | — |
| `calcHourglassResist.f90` | FOLDED — `assembleGlobalKU.py:471-515`, a separate function called separately from `driver.py:140` | — | — |
| `calcLocalShapeFunc.f90` | FOLDED — `assembleGlobalMass.py:111-124` (`_LOCAL_DERIV`, `_W`) | — | — |
| `calcQAttenuationCoeff.f90` | **ABSENT** | no Q attenuation; no `C_Q` in `readInputFiles.py:43-128` | **UNREACHABLE ON BOTH SIDES**: `C_Q` is hardcoded `0` at `globalvar.f90:104` and read from no input file, so `calcElemKU.f90:77`'s branch cannot execute without editing source. Not a port gap |
| `checkInputConsistency.f90` | `checkInputConsistency.py` (`check`), called from `eqdyna3d.py`'s `build_solver_state` right after `readInputFiles.build_params` — the same point `eqdyna3d.f90:104` calls it, before any mesh or solver work. Raises `InputConsistencyError` carrying the matching `errorCodes.f90` `ERR_CFG_*` number; `eqdyna3d.py main()`'s `_abort` turns that into `SystemExit(code)` on both the serial and `--mpi` paths, so `python3 -m eqdyna` exits with the same code a Fortran run of the same bad config exits with. Two of the three ported checks (`C_Q==1` combinations) stay unreachable through any real case on BOTH sides — `C_Q` is never case input (see the `calcQAttenuationCoeff.f90` row) — and are pinned by direct unit call in `testsys/regression/test_check_input_consistency.py`, not through a case directory | — | — |
| `computePMLDampingVector.f90` | FOLDED — `func_lib.py:59-62` (the 3·vmax/2δ·log(1/R)·(d/δ)² scaling, which Fortran keeps OUT of `pmlRegionDistance`) + `assembleGlobalKU.py:133-137` (the `dv(1:9)` tiling). The negative-component guard (`:40-44`) is not folded | — | — |
| `countMeshEntities.f90` | FOLDED — `meshgen.py:323-408` (nodes, `nftnd`), `:410-646` (elements incl. `wedge4num`'s extra count), `:898-1026` (equations). Counts are array shapes; one pass instead of two | — | — |
| `driver.f90` | `driver.py` | — | — |
| `eqdyna3d.f90` | `eqdyna3d.py` | — | — |
| `errorCodes.f90` | **PARTIAL as a numbered-code contract** — the `ERR_CFG_*` (11-13) range now has a Python analogue: `checkInputConsistency.py`'s `InputConsistencyError.code` + `eqdyna3d.py`'s `_abort`, matched code-for-code against this registry. Every other range (`ERR_INPUT_FILE_*`, `ERR_GEOM_*`, `ERR_MESH_*`, `ERR_MPI_*`, `ERR_NUM_*`, `ERR_NETCDF`) still has no Python counterpart | outside `ERR_CFG_*`: `ValueError`/`NotImplementedError`, no numbered exit codes and no `abortRun` status contract | NOTHING GATED for the exit-code contract outside `ERR_CFG_*`. **But see the re-verification below: the Jacobian guard the audit reported as missing DOES exist in the port.** `ERR_NUM_PML_DAMPING` and `ERR_NUM_VELOCITY_NAN` analogues are unverified either way |
| `faulting.f90` | `faulting.py` | — | — |
| `fric.f90` | `fric.py` | — | — |
| `func_lib.f90` | `func_lib.py` | — | — |
| `globalvar.f90` | `globalvar.py` | — | — |
| `library.f90` | **ABSENT** (`memory_estimate`) | the port prints no memory estimate | NOTHING GATED — stdout is not compared, the sweep compares `frt` only. Zero physics |
| `library_degeneration.f90` | FOLDED — `meshgen.py:174-185`, `:186-206`, `:540-578`, `:648-664`, `:666-725`, `:831-870`. **`C_degen>3` IS ported**; `test.tpv36`/`test.tpv37`'s python-numpy and python-jax cells run the real wedge path, gated at 1e-6, and are in `CI_CELLS` | — | — |
| `library_output.f90` | `library_output.py` | — | — |
| `meshgen.f90` | `meshgen.py` | — | — |
| `netcdf_io.f90` | **PARTIAL** — `netcdf_read_on_fault_eqdyna` FOLDED into `readInputFiles.py:300-347`; `netcdf_read_on_fault_eqdyna_restart` (`:114-185`) **ABSENT** | restart `*.r.nc` is never read; the port always initialises from `on_fault_vars_input.nc` | NOTHING GATED — no gated case restarts from a previous cycle |
| `readInputFiles.f90` | `readInputFiles.py` | — | — |
| `updateThermalPressurization.f90` | `updateThermalPressurization.py` | — | — |

## `CLAUDE.md` drift found by the audit

- **D1** — CLAUDE.md says `assembleGlobalKU.py` covers four Fortran files
  "because they are one fused loop in the port". `driver.py:139-140` calls
  `KU.assembleGlobalKU(...)` and `KU.calcHourglassResist(...)` as two separately
  defined functions (`assembleGlobalKU.py:261`, `:471`). The mapping is right;
  the stated reason is wrong — they are co-located, not fused.
- **D2** — the same sentence claims `calcElemMass` whole; the module's own
  docstring says "`calcElemMass.f90`'s gravity body force", i.e. partial. The
  `rdampm` half is not carried. **Resolved by the conductor, see below.**
- **D3** — "every Python module is named after its Fortran counterpart" with
  four stated exceptions undercounts: `meshgen.py` also carries
  `library_degeneration` and `countMeshEntities`, `readInputFiles.py` also
  carries `netcdf_io`'s reader, `func_lib.py` also carries
  `computePMLDampingVector`'s damping-profile scaling. Seven exceptions, not
  four. This table is the fix.
- **D4** — `func_lib.py:4-6` says `region_damp` "returns the three damping
  coefficients"; `func_lib.f90:16-77` returns raw DISTANCES and both Fortran
  callers apply the scaling themselves. Comment-vs-code drift, and the reason
  `computePMLDampingVector` looks absent when it is folded.

## Conductor's re-verification (wei-lin, 2026-09-23) — one finding refuted, one resolved

**REFUTED: "an inverted element makes Fortran stop and Python produce a finite,
wrong answer."** The audit's top-ranked absence does not hold. Both of the
port's determinant paths refuse a non-positive Jacobian:

- `assembleGlobalMass.py:147-167`, `compute_element_det` — docstring states it
  outright ("Raises ValueError on any non-positive determinant (Fortran's own
  `if (det <= 0.0d0) ... stop`, made a loud Python...)"), and `:167` is
  `bad = np.nonzero(det <= 0.0)[0]`.
- `assembleGlobalMass.py:461-464`, `compute_element_shape` — the path actually
  used by `eqdyna3d.py:301` — raises `ValueError('compute_element_shape:
  non-positive determinant at element(s) ...')`.

What survives of that row is narrower and still true: the port has no NUMBERED
exit-code contract, so a refusal is a Python traceback rather than a code from
`errorCodes.f90`'s registry. That is a real declared absence; "no Jacobian
check" is not.

**RESOLVED: D2 is a documentation fix, not a physics gap.** `rdampm` appears
exactly twice in `src/fortran/` — the declaration `rdampm=0.0d0`
(`globalvar.f90:179`) and its single use (`assembleGlobalKU.f90:19`). No input
file and no reader sets it, so the Fortran value is 0 for every case that can
be run, and the port's hardcoded `rdampm=0.0` (`eqdyna3d.py:358`) matches the
Fortran exactly. Mass-proportional Rayleigh damping is unreachable on both
sides, like `C_Q`.

**Standing, unverified by me** (recorded as the audit reported them, not
independently re-derived): the `computePMLDampingVector` negative-damping guard
and the velocity-NaN guard having no port analogue; the restart-read absence.
Each is a candidate board row. (`checkInputConsistency`'s absence, also listed
here by the original audit, is CLOSED — see its row above.)

**One folding stands on no measurement.** `calcB` — never forming B — is
nowhere justified, measured or otherwise. Under rule 23 clause 5 that makes it
a departure to be recorded under its real reason, not a perf claim. The other
foldings do carry measurements in the source (`assembleGlobalMass.py:57-71`,
`assembleGlobalKU.py:274-286` and `:14-27`).

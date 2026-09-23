# Item 71 — coverage for the generated (`insertFaultType=2`) fault surface

Worktree `/home/utig5/dliu/eqdyna-wt-item71`, branch `item71-fractal-coverage`,
from `5b7a278` (v5.16.1). 2026-09-23, iris-vermeulen.

## What was added

`testsys/regression/test_fractal_fault_geometry_derivatives.py` — a SIBLING of
`test_rough_fault_normal_consistency.py`, not an extension of it.

Why a sibling: that file is the incident record for the SHIPPED SCEC surface
(`702b56b`). Its fixture (`tpv29_official_25m_block.txt`), its constants
(`OFFICIAL_DX`, `GATE_DX`, TPV29/30's R / a / B11 / B33 / B13) and all six of
its checks are TPV29/TPV30's, and it reaches them through `case_input/` on
disk. The generated path shares none of that: it needs `matplotlib`, a
module-level `par` stub, and a tmpdir. Two files means a failing banner names
the path that broke. `testsys/run.py:run_regression` globs `test_*.py`, so the
sibling is registered by existing.

## The assertions

1. `generated_derivatives_at_own_spacing` — the one item 71 asked for. A fresh
   `insertFaultType=2` surface (32x16, dx=100, dz=50, dip=90, seedId=1 fixed by
   the script) has stored `dy/dx` == `np.gradient(y, dx, axis=1)` and stored
   `dy/dz` == `np.gradient(y, dz, axis=0)`, max|diff| 5.000e-07 against a
   1.0e-4 absolute bound. `dx != dz` on purpose: a spacing swap cannot pass.
2. `fractal_branch_produces_roughness` — a plane satisfies (1) trivially, so
   (1) alone would pass against a dead type-2 branch. Type 2 gives max|slope|
   0.4505; type 1 at dip 90 gives 0.0 and a flat surface.
3. `wrong_spacing_is_discriminated` — differentiating the same surface at 2dx /
   2dz misses by 0.2253 = 2253x the bound.
4. `production_gate_accepts_generated_surface` — `lib.
   validateFaultRoughGeometryForCase`, the exact call at `case.setup:379-380`,
   run on the generated file WITHOUT case.setup.
5. `identity_refuses_perturbed_derivative` — +0.05 on one stored `dy/dx` puts
   check 1's residual at 510x its bound.
6. `production_gate_refuses_perturbed_derivative` — +0.2 on the same value is
   refused by the production gate.

## Measurement worth keeping (it set the shape of checks 5/6)

Bisection on this surface, 2026-09-23: the smallest single-node `dy/dx` error
`validateFaultRoughGeometryForCase` refuses at the probed node is **0.0568** —
12.6% of the surface's own max|slope| of 0.4505. A +0.05 perturbation is
ACCEPTED by it. That is `lib._checkDerivativeColumn`'s data-implied truncation
bound behaving exactly as `lib.py:556` documents (a SUPPLIED file may carry the
analytic derivative, so the tolerance must be derived from the data), on a
surface deliberately near the generator's own 0.2 roughness warning. It is not
a defect in `lib`. It IS the argument for this guard: a GENERATED file writes
`np.gradient` itself (`generateFaultInterface:94-95`), so it can be held to an
absolute 1e-4 — ~570x tighter than the gate that backs it up.

## Mutation test (run, reverted, not committed)

Temporary two lines after generation moved ONE stored `dy/dx` value by 0.001 —
deliberately UNDER the 0.0568 gate threshold:

```
generated_derivatives_at_own_spacing: FAIL - ... max|diff| 9.995e-04 ... bound 1.0e-04
production_gate_accepts_generated_surface: SUCCESS - ...
fractal_fault_geometry_derivatives: FAIL - 1 check(s) failed: generated_derivatives_at_own_spacing
```

The new check catches a stored-column regression that `case.setup`'s own gate,
on the same file, still accepts. Reverted; the file re-runs at exit 0.

## Not done / for the row's author

- `pathway_forward.md` item 71's verification Command is
  `grep -c insertFaultType testsys/regression/test_rough_fault_normal_consistency.py`.
  It is hardwired to ONE file and still reads 0 with this guard landed, i.e. it
  reports "still uncovered" when it is covered. It should be
  `grep -rl insertFaultType testsys/regression/` (reads the new file) or name
  `testsys/regression/test_fractal_fault_geometry_derivatives.py` directly.
- No frt.canonical.txt, no matrix.py bound, no solver source, no
  PROJECT_RULES.md, no pathway_forward.md touched.

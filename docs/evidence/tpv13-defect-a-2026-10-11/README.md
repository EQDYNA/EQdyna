# TPV13 defect A: before / after (2026-10-11)

Owner decision 2026-10-11 (board PR #178): the next small release carries the
fix for TPV13 defect A, named in `docs/evidence/barall-overlays-2026-10-09/resolution/`.

## Root cause

`createElement` (`src/fortran/meshgen.f90`, port: `src/python/eqdyna/meshgen.py`
`build_elements`) typed each element PML or regular by the mean of its
**y-blended** `meshCoor`. `setNumDof` gives a node 12 dofs from the
**unblended** grid lines. On dipping cases (tpv10/12/13) `insertFaultInterface`'s
y-blend shifts the footwall by as much as `peak` = 7.5 km toward the fault, so
near the y-min PML interface regular (type-1) elements ended up holding 12-dof
nodes. Their y/z element force went into the split PML x-dofs, so the Method-2
(`C_elastic=0`) initial stress at the y-min PML interface was out of balance.
The result was a spurious P wave launched at t=0 from the y-min face.

Fix: type an element by the centre of its grid-line cell (`gridCenterCoor`;
the port uses `0.5*(xl[IX-1]+xl[IX])` etc.). Changed in both implementations.

## Measurement that proves it

The nt==1 equilibrium residual (`EQDYNA_DUMP_EQUIL=1`; `dumpNodalAccel` in
`driver.f90` now sums a 12-dof node's split dofs) for test.tpv13 at dx 500 m:

| | max nodal accel at nt==1 |
|---|---|
| before | ~90 m/s^2, one slanted node column at the y-min PML interface |
| after | 2.944e-2 m/s^2, the generic free-surface value |

Census on a coarse tpv13 mesh (dx 1500 m, 74,400 elements, 42,450 12-dof
nodes): regular elements holding a 12-dof node, 798 before and 0 after.
`testsys/regression/test_pml_element_type_grid_centre.py` checks this and
comes out both ways: SUCCESS on this branch, FAIL on `fd56af1`.

## Outcome (falsified against FaultMod, not only against the defect)

These are `make_overlays.py resolution` figures and tables (the
`scripts/figures/scec_compare.py` overlay engine) at dx 500 m against FaultMod.
`before/` is `fd56af1`; `after/` is this fix.

| tpv13, 8 s window | before | after |
|---|---|---|
| ruptured | 91.2% | 94.2% |
| p90 \|dt\| (s) | 2.838 | 0.531 |
| max \|dt\| (s) | 3.72 | 2.67 |
| slip ratio | 0.61 | 0.95 |
| body peak-velocity ratio | 4.24 | 1.01 |

`*_ts_body_res.png`: before, EQdyna leaves FaultMod after 3-4 s, with
v-velocity down to -10 m/s and n-displacement +25-32 m against FaultMod's
~2-3 m/s and ~-6 m. After, the curves track FaultMod. The y-min boundary wave
is gone; the residual table above is its source term.

## TPV27

Not affected. On master the nt==1 residual is <= 4e-12 m/s^2, because the
fault is vertical and nothing is blended toward a PML face. The Method-2
free-surface sinking was defect C, fixed in PR #167.

## Which references move

Defect A meshes exist in test.tpv10, test.tpv12 and test.tpv13, which share the
dipping geometry. But tpv10 and tpv12 are Method 1 (`C_elastic=1`, no initial
body stress), so the mistyped elements carry no unbalanced stress. At the 5 s
gate term both still pass their committed references on both backends: worst
station e is 5.1e-09 (tpv10) and 4.5e-14 (tpv12). Only test.tpv13 (Method 2)
moves: worst station e is 1.62 against bound 1e-10. Its references are
regenerated in their own commit (rule 7).

## Not in this change: the TPV13 border unpin and the fault-edge master node

The SCEC TPV13 spec (p.3) puts no strength barrier at the fault border, so the
mu_s=1000 pin in `case_input/test.tpv13` should go. Removing the pin on top of
this fix exposed a second, independent defect. `replaceSlaveWithMasterNode`'s
+x test `nodeCoor(1) < fxmax+dx+tol` assumes the cell past `fxmax` is exactly
dx wide. On a graded grid (15000 -> 15512.5 m) the hanging-wall element past
+fxmax keeps the slave edge node. Edge tractions become x-asymmetric, and in
Method 2 the +x edge breaks acausally at t=0 (|tau|/|sigma_n| 0.81 > mu_s 0.70).
The fix tests the cell's lower grid line (`xline(ix-1) < fxmax+tol`), in both
implementations, on branch `fix/fault-edge-master-node`. With defect A, the
edge fix and the unpin together: 97.1% ruptured, 0 acausal nodes, slip ratio
0.96, body pv ratio 1.01, 58/121 border nodes ruptured (earliest 4.09 s).
That fix moves the references of nearly every registered case, so it follows
as its own PR with a suite-wide regeneration.

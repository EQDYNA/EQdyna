# TPV13 fine-grid (100 m) FaultMod overlay, 2026-10-11

Owner acceptance criterion: a released TPV must have a fine-grid run (typically 100 m)
that matches its published / Barall independent result on the same grid, judged by
visual cplot + ts overlay. This is that run for test.tpv13.

## Run

- `test.tpv13`, `par.dx = 100` m (spec resolution), `par.term = 8` s (spec term,
  `testsys/e2e/full_specs.py` `test.tpv13` entry), EQdyna fortran, 8 MPI ranks
  (`par.nx,ny,nz = 4,2,1`; bumped from the committed `2,2,1` for wall-clock only --
  this is a scratch run, no committed reference touched, rule 7). Wall time: 80 min.
- Baseline: `~/shared_dataset/scec_cvws.tpv1213/raw/tpv13/barall-faultmod-100m-2009`
  (Barall FaultMod, also 100 m) -- the SAME grid, so this is an exact-node comparison
  with nothing interpolated (45,280 shared fault nodes at 8 s, vs 1,877 at the old
  gate dx 500 m where only a coarse subset of FaultMod's 100 m grid overlapped).
- Built with the existing `docs/evidence/barall-overlays-2026-10-09/make_overlays.py
  resolution` tool and `scripts/figures/scec_compare.py`; only new input is
  `RUNS/tpv13_dx100/`.
- Code state: master @ d225da7, i.e. **Defect A fixed** (PR #180, grid-line-centre PML
  typing) and **Defect C fixed** (PR #167, +7.3215 m depth offset removed). **Defect B
  (fault-border pin, item 155) is NOT fixed** -- `fix/fault-edge-master-node` (d71f3b0)
  stays unmerged per owner decision 2b (2026-10-11); this run does not touch it.

## Verdict: MATCH by the declared rupture/slip criteria (`../barall-overlays-2026-10-09/README.md`); NOT yet accepted under the owner's visual criterion

The owner criterion (board, PR #185) also requires the time series to visibly match.
They do at every buried station and for every slip / slip-rate / shear-stress column;
the free-surface fault node's normal stress (and, through friction, its v-shear
stress) does not -- see "Surface-node normal stress" below. That mismatch is real,
measured, and not yet explained; TPV13 stays open under the owner criterion.

| window | dx | ruptured | median \|dt\| | p90 | one-code-only | slip ratio | body pv ratio | MATCH |
|---|---|---|---|---|---|---|---|---|
| 5 s | 100 m | 99.6% (43575/43730) | 0.032 s | 0.067 s | 0 | 1.01 | 1.02 | yes |
| 8 s (spec) | 100 m | 98.8% (44756/45280) | 0.034 s | 0.070 s | 0 | 0.99 | 1.02 | yes |

Criteria: ruptured >= 95%, median <= 0.2 s, p90 <= 0.5 s, no one-code-only station --
all four hold at both windows, by a wide margin (median/p90 are ~6-7x inside bound).

This resolves the open "DEFECT" verdict against this same baseline recorded in
`../barall-overlays-2026-10-09/resolution/README.md` (there: 125 m, pre-Defect-A/C,
slip ratio 0.79, body pv ratio 3.58, "yes (rupture timing only)"). That table is now
superseded for the fine-grid question: it predates both fixes (committed 2026-10-09,
before PR #167 2026-10-10 and PR #180 2026-10-11), and this is the first 100 m run on
top of both.

## Visual read (PNGs, read directly, not just the numbers)

- `tpv13_cplot_res.png`: rupture-front contours overlap almost exactly across the
  whole fault, 1 s intervals, 0-8 s. The only visible separation is at the fault
  BORDER (|x| near 15 km and depth near 13-15 km), where ours lags / stays inside
  FaultMod's contour -- the expected signature of the still-pinned Defect B border
  (`mu_s = 1000` at |x|=15 km and the 15 km down-dip row), not a new defect.
- `tpv13_ts_fault_res.png` (regenerated 2026-10-12 with the spec station names, see
  "Station naming" below; 10 of 10 discovered on-fault stations now comparable, was 3
  of 17 with 14 excluded by name): h-slip, h-slip-rate, h-shear-stress, v-slip,
  v-slip-rate overlay near pixel-for-pixel at every station shown, and so do
  v-shear-stress and n-normal-stress at every BURIED station (dp015 -24.2 vs FaultMod
  -25.4 MPa at 8 s, dp030 -36.9 vs -36.7, dp045 -42.3 vs -42.5, dp075 -58.3 vs -58.6,
  dp120 -88.6 vs -88.5). **At the three free-surface stations (dp000) n-normal-stress
  does not overlay**, and v-shear-stress follows it (below).
- `tpv13_ts_body_res.png`: all four off-fault stations shown (h/v/n displacement and
  velocity) overlay closely; the largest departure is v-displacement at
  body-030st000dp000, within a few percent.

## Not done in this pass

- SCEC CVWS submission zip: no existing builder script in this repo
  (`case_input/test.tpv13/README.md` references the comparison test, not a packager).
  Building one correctly needs the SCEC submission format spec read in full, which
  this pass did not budget for. Left as a follow-up, not attempted partially.
- The free-surface n-normal-stress divergence: investigated (below), cause not found;
  routed as a physical-validity question (corner traction at the fault / free-surface
  intersection of a 60-degree normal fault), not patched with an invented
  regularisation.

## Station naming (fixed, PR #188)

`library_output.f90` divided the station depth by `dsin(fltxyz(2,4,j))`, which is
forced to 90 degrees for `insertFaultType > 0, C_degen == 0`, so tpv10/12/13 stations
were named by vertical DEPTH (`dp013/dp026/dp039/dp065/dp104`) rather than the spec's
down-dip distance (`dp015/dp030/dp045/dp075/dp120`). The reader's name heuristic
rescued the 500 m files but not the 100 m ones (14 of 17 excluded in the first version
of these overlays). PR #188 fixes Fortran + port and renames the gated references;
the overlays here were regenerated from the same run with the files renamed (solver
output bytes unchanged, names and `# location` header only).

## Surface-node normal stress (investigated 2026-10-12, open)

n-normal stress, MPa, sign per spec output (negative = compression), at 4 s / 8 s,
from `faultst*dp000.txt` of this run and FaultMod's:

| station | EQdyna 100 m | FaultMod 100 m | EQdyna 200 m | EQdyna 500 m |
|---|---|---|---|---|
| faultst000dp000 | +8.43 / +10.01 | +1.47 / +2.93 | +1.96 / +3.31 | -- / -11.74 |
| faultst045dp000 | +3.04 / +3.77 | -0.71 / -0.87 | -- | -- / -12.34 |
| faultst120dp000 | -0.17 / -0.17 | -0.16 / -0.91 | -1.78 / -5.15 | -- / -14.54 |

(200 m and 500 m: 8 s, 4 ranks, same code; 500 m is the gate grid.)

- v-shear stress is NOT an independent failure: with `mu_d = 0.1`, cohesion 0.2 MPa and
  tension clamped to zero effective normal stress, post-rupture shear = 0.2 + 0.1 x
  max(-sigma_n, 0). FaultMod st045 at 8 s: 0.287 = 0.2 + 0.1 x 0.87; ours: 0.200
  because our surface node is tensile. st000: both 0.2000. Fix the normal stress and
  the shear column follows.
- Not a sign/plot convention: both codes report -0.369 MPa at t=0 at the surface
  station (7390.01 Pa/m x 50 m down-dip), and the buried stations agree to <5%.
- Not hourglass control or Rayleigh damping: at 500 m / 8 s, baseline -11.74 /
  -12.34 / -14.54 (st000 / st045 / st120); viscous hourglass (`C_hg = 2`) -9.88 /
  -10.48 / -10.86; `rdampk = 0` -11.41 / -12.00 / -15.26 (<20% change, wrong direction
  to reach FaultMod).
- Confined to the free-surface row: down-dip profile at x = 0, 8 s, 100 m: +10.0,
  -5.0, ... converging to the buried values within ~1.5 km; element stresses adjacent
  to that node (`pstr.txt`) are small and compressive (sigma_yy -1.87 footwall / -0.57
  hanging wall MPa), so the +10 MPa is a nodal-traction quantity not reflected in the
  element field.
- Mesh-non-convergent in ours (500 m -11.7, 200 m +3.3, 100 m +10.0, sign-changing,
  smooth bell along strike at 100 m, +10 at x = 0 to ~0 at |x| = 8 km), while the
  rows below converge. Interpretation: the fault / free-surface corner of a dipping
  frictional interface carries a discretisation-dependent traction in the TSN
  split-node formulation; FaultMod's +2.9 is its own discretisation's value. Whether
  either is "right" is a physical-validity question, not a code-parity one.
- Run inputs for the comparison: `RUNS/tpv13_dx100` (this run), `barall-faultmod-100m-2009`
  (read-only baseline); the 200 / 500 m runs and the hourglass / damping spikes are
  scratch (`.claude/worktrees/tpv13-surface/RUNS/`, not committed).

## Reproduce

```
export EQDYNAROOT=<worktree>; export PATH=$EQDYNAROOT/bin:$EQDYNAROOT/scripts:$PATH
scripts/create.newcase RUNS/tpv13_dx100 test.tpv13
# edit RUNS/tpv13_dx100/user_defined_params.py: par.dx = 100.
(cd RUNS/tpv13_dx100 && ./case.setup && mpirun -np 8 eqdyna)
python3 -m testsys.frt_canonical RUNS/tpv13_dx100
python3 docs/evidence/barall-overlays-2026-10-09/make_overlays.py resolution RUNS \
    docs/evidence/fine-grid-overlays-2026-10-11/tpv13 tpv13
```

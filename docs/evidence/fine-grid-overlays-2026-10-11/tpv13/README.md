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

## Verdict: MATCH (both windows), by the declared criteria (`../barall-overlays-2026-10-09/README.md`)

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
- `tpv13_ts_fault_res.png`: h-slip, h-slip-rate, h-shear-stress, v-slip, v-slip-rate
  and v-shear-stress all overlay near pixel-for-pixel at all three stations shown.
  **n-normal-stress does not**: at faultst000dp000 and faultst045dp000 ours settles
  near +10 MPa post-rupture vs FaultMod's ~0-3 MPa. Initial (t=0) normal stress
  matches FaultMod to 5 significant figures (checked separately in the defect-A
  evidence's station table), so this is a POST-RUPTURE divergence, not a setup bug.
  It is not scored by the declared MATCH criteria (which are slip/rate/timing-based),
  and the suite's own gate design (CLAUDE.md "There is ONE test") checks n-stress by
  SIGN at buried stations, not magnitude -- this is flagged here as an open question
  for whoever looks at TPV13 plasticity next, not chased in this pass.
- `tpv13_ts_body_res.png`: all four off-fault stations shown (h/v/n displacement and
  velocity) overlay closely; the largest departure is v-displacement at
  body-030st000dp000, within a few percent.

## Not done in this pass

- SCEC CVWS submission zip: no existing builder script in this repo
  (`case_input/test.tpv13/README.md` references the comparison test, not a packager).
  Building one correctly needs the SCEC submission format spec read in full, which
  this pass did not budget for. Left as a follow-up, not attempted partially.
- The n-normal-stress post-rupture divergence noted above: measured, not diagnosed.

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

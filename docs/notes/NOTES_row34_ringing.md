# Row 34 -- post-front on-fault stress ringing at 500 m (tpv36/tpv37), controlled experiment 2026-09-29

Owner's target (2026-09-29): `faultst000dp180`, down-dip shear stress
(column 7). The 500 m gate reference reaches -34.3 MPa by ~0.5 s and then
zigzags; the published 50 m SCEC run
(`scec_archive/tpv37/eqdyna-v5.3.3-50m-2024-nstress-corrected/faultst000dp180`,
read-only) is flat at the same level.

**Verdict: no hourglass/damping setting removes the zigzag without moving
the rupture. The zigzag is a resolution artefact: it falls by 3x per halving
of dx (rms 0.150 -> 0.051 MPa at 500 -> 250 m) and the 50 m archive sits at
0.0018 MPa on the same metric. The one non-destructive damping change found
(KF78 + viscous, kapa_hg = 0.03) buys 45% in rms and 21% in peak-to-peak
while shifting 18/3477 rupture flags and peak slip rate by -11%: not a fix.**

## Metric (`testsys/perf/measure_ringing.py`, committed with these notes)

Per on-fault station, on `v-shear-stress`:

- `t_arr`: first sample with `hypot(h-slip-rate, v-slip-rate) > 1e-3 m/s`
  (SCEC convention). Never reached -> NOT-RUPTURED, does not vote.
- window `[t_arr + 0.5 s, min(t_arr + 2.5 s, t_end)]`; refused if < 1.0 s.
  The 0.5 s lead excludes the stress drop (grid period at 500 m is 0.29 s;
  the drop at these stations completes within ~0.3 s of arrival).
- cubic least-squares detrend over the window; ringing = RMS and
  peak-to-peak (MPa) of the residual, plus its dominant frequency.
- `max|run - other|` after `np.interp` of `other` onto the run's time axis.
- Rupture: run `frt.txt*` canonicalised and aligned to the committed
  `frt.canonical.txt`: ruptured-node counts, flips, `|d rupture time|`
  max/mean, fault-wide peak slip rate (frt col 10) ratio.
- `--max-p2p X`: exit 1 (UNFIXED) if any measurable station exceeds X, 0
  (FIXED) otherwise, 2 if no station is measurable. Rule 14a: the committed
  reference reads UNFIXED at 0.2 MPa today (p2p 0.673); the 50 m archive
  file, fed as `--station`, reads FIXED at the same bound (p2p 0.016).

At the 5 s gate term two of tpv37's three on-fault GATE stations
(`faultst000dp010`, `faultst120dp030`) have not ruptured, so the table
below uses dp180 plus `faultst000dp120` / `faultst040dp180`, which do.

## Runs

All Fortran, 4 ranks (`taskset -c 2-5|6-9 mpirun --bind-to none -np 4`),
`par.term = 5.0`, dx = 500 m unless stated, worktree `research/ringing-36-37`
at `c9531d6`, cotopaxi (AMD EPYC 7532, load average 85-100 throughout:
timings are not reported, physics is deterministic), gfortran 11.4.0,
Open MPI 4.1.1, `OMP_NUM_THREADS=1`. `C_hg`/`kapa_hg` are compiled in
(`src/fortran/globalvar.f90:105,191`); each variant is a scratch copy of
`src/fortran` with those two lines changed, built by the same makefile.
`C_hg=3` below is a scratch-only edit of `calcHourglassResist.f90` making
BOTH branches execute (KF78 stiffness + Goudreau-Hallquist viscous); it
does not exist in the tree. `rdampk` is `par.rdampk` (default 0.1, x dt at
read time). Nothing under `src/` in the worktree was changed.

Sanity: `tpv37.default` and `tpv36.default` reproduce their committed
references at max|diff| = 0.0 over all 22 frt columns, and the station
series match to 0.0000 MPa.

### tpv37, `faultst000dp180`, down-dip shear stress, window [0.81, 2.81] s

| config | ringing rms / p2p (MPa) | f_dom (Hz) | max\|run - 50 m archive\| (MPa) | ruptured nodes (ref 1288 / 3477) | flips | \|d rupt. time\| mean / max (s) | fault peak slip rate ratio | dp180 slip at 5 s (m, ref 1.213) |
|---|---|---|---|---|---|---|---|---|
| **default (C_hg=1, rdampk=0.1)** | 0.150 / 0.673 | 7.5 (broadband 4.5-8) | 1.246 | 1288 | 0 | 0 / 0 | 1.000 | 1.213 |
| 50 m archive (v5.3.3) | 0.0018 / 0.016 | -- | 0 | -- | -- | -- | -- | 3.272 |
| dx = 250 m, default | 0.051 / 0.402 | 16 | 0.743 | not comparable (different mesh) | -- | -- | -- | 1.324 |
| C_hg=2, kapa_hg=0 (no HG control) | 0.245 / 1.585 | 1.0 | 2.150 | 614 | 674 | 0.364 / 2.51 | 1.305 | 0.920 |
| C_hg=2, kapa_hg=0.01 | 0.203 / 1.387 | 1.0 | 1.875 | 541 | 747 | 0.439 / 2.52 | 1.153 | 0.857 |
| C_hg=2, kapa_hg=0.03 | 0.143 / 1.047 | 0.5 | 1.393 | 520 | 768 | 0.559 / 2.79 | 0.943 | 0.799 |
| C_hg=2, kapa_hg=0.1 (default coeff.) | 0.039 / 0.302 | 0.5 | 0.537 | 436 | 852 | 0.560 / 2.47 | 0.565 | 0.745 |
| C_hg=3 (KF78 + viscous), kapa_hg=0.03 | 0.083 / 0.529 | 5.0 | 1.253 | 1270 | 18 | 0.015 / 0.20 | 0.886 | 1.204 |
| C_hg=3, kapa_hg=0.1 | unstable (55 MPa p2p, peak slip rate 5.6e4 m/s) | | | 3477 | 2189 | | 11149 | |
| rdampk = 0.3 | unstable (values ~1e22 by 5 s) | | | 3477 | 2189 | | | |
| rdampk = 1.0 | aborted: "Velocity became NaN during time stepping" | | | -- | | | | |

Other stations, default vs dx=250 (rms / p2p MPa): `faultst000dp120`
0.061 / 0.247 -> 0.017 / 0.097 (archive 0.024 / 0.123 -- at this station
the 250 m run is already at the archive's level); `faultst040dp180`
0.112 / 0.519 -> 0.061 / 0.311 (archive 0.008 / 0.046).

tpv36 default at dp180 is numerically identical to tpv37's at 5 s
(0.150 / 0.673; the two cases differ only in cohesion slope and the
hypocentre station has not felt it by then). Not re-run for the variants.

### Control on a hex mesh: tpv10, C_hg=2, kapa_hg=0.1

`faultst000dp104` ringing 1.624 / 8.611 -> 0.939 / 3.465 MPa, but ruptured
nodes 714 -> 141 (573 flips of 1891), peak slip rate ratio 0.962. Same
shape as tpv37: the rupture collapse under C_hg=2 is not a wedge
(C_degen) effect.

## What the sweep says

1. **Switching to `C_hg=2` collapses the rupture almost independently of
   `kapa_hg`**: ruptured nodes 614 / 541 / 520 / 436 at kapa 0 / 0.01 /
   0.03 / 0.1. The kapa = 0 run has NO hourglass control at all and already
   loses 52% of the ruptured area. So the board's reading -- "C_hg=2
   over-damps the rupture" -- is mostly wrong in mechanism: the bulk of the
   damage is the LOSS of the KF78 stiffness control (`elseif` makes the two
   branches exclusive), and the viscous term adds to it (614 -> 436). The
   viscous term alone also does not hold the hourglass modes: ringing is
   WORSE than default at kapa <= 0.03 (1.05-1.59 vs 0.67 MPa p2p).
2. **KF78 + viscous (C_hg=3) is the only non-destructive damping found**, at
   kapa_hg = 0.03: -45% rms, -21% p2p, 18 flips, mean rupture-time shift
   0.015 s, peak slip rate -11%, dp180 slip -0.7%. Still 46x the archive's
   rms. At kapa_hg = 0.1 it is unstable at the case's CFL of 0.5
   (`par.dt = 0.5*dz/vp`, and dz = dx sin 15 deg = 129 m is the stiff
   direction).
3. **`rdampk` cannot be raised**: 0.3 blows up, 1.0 NaNs. The default 0.1
   is close to the explicit scheme's stability limit for this anisotropic
   mesh, so the "0.65% of critical" damping the board computed cannot be
   increased by this knob.
4. **Refinement does what the board predicted it would not do
   uniformly**: rms 0.150 -> 0.051 at 500 -> 250 m (0.34x per halving),
   and 0.0018 at 50 m; extrapolating 0.34x over the remaining 2.3 halvings
   gives ~0.004, within 2x of the archive. The frequency moves up with
   Vs/(2dx) (7.5 -> 16 Hz). The board's "resolution-independent zeta"
   argument is about the DECAY RATE of the mode; the AMPLITUDE the rupture
   front puts into it falls with dx. The mesh is the actor.
5. Error vs the 50 m archive at dp180 is dominated by the front (arrival /
   drop shape), not the ringing: default 1.25, C_hg=3 k0.03 1.25 (ringing
   halved, error unchanged), dx=250 0.74.

## Tried and rejected

| Approach | Why it failed | Evidence |
|---|---|---|
| C_hg=2 at its default kapa_hg=0.1 (board's lead) | ringing -55% but rupture area -66%, peak slip rate -43% | table row; tpv10 control same shape |
| C_hg=2 with smaller kapa_hg (0.03, 0.01) | rupture area still -60%, ringing WORSE than default | 520 / 541 nodes vs 1288; 1.05 / 1.39 MPa p2p |
| C_hg=2 with kapa_hg=0 (mechanism control) | isolates the KF78 loss as the main actor | 614 nodes, 1.59 MPa p2p |
| KF78 + viscous, kapa_hg=0.1 | unstable at CFL 0.5 | 55 MPa p2p, 5.6e4 m/s |
| KF78 + viscous, kapa_hg=0.03 | works but buys 45% rms, not the 50x needed; -11% peak slip rate | table row |
| rdampk 0.3 / 1.0 | explicit-scheme instability | 1e22 values / NaN abort |

## What a fix would touch, if one were wanted anyway

A `C_hg=3` (stiffness + viscous) with a runtime `kapa_hg` would touch
`src/fortran/globalvar.f90` (both scalars are compiled in today),
`readInputFiles.f90` + `scripts/case.setup` + `scripts/defaultParameters.py`
(to make them case parameters), `calcHourglassResist.f90` (branch
structure), and the port at `src/python/eqdyna/assembleGlobalKU.py`
(`calcHourglassResist` there is KF78-only; a viscous term is a new scatter
contribution in the fused loop, rule 23 state (b) record to update in
`docs/fortran_python_correspondence.md`). Not recommended on this evidence:
the only stable, non-destructive coefficient leaves the zigzag at 46x the
archive, and it changes every committed 500 m reference (rule 7 commit).

The honest route to a flat dp180 trace is a finer gate mesh, and that is a
suite-cost decision (dx=250 tpv37 at 4 ranks took ~28 min on this loaded
box vs ~3 min at 500 m, i.e. ~9x, elements x steps = 16x nominal).

## Reproduce

```
# baseline / any run dir against reference + archive
python3 testsys/perf/measure_ringing.py --case test.tpv37 --run <run_dir> \
  --archive scec_archive/tpv37/eqdyna-v5.3.3-50m-2024-nstress-corrected \
  --stations faultst000dp180,faultst000dp120,faultst040dp180 --max-p2p 0.2
# committed reference alone: reads UNFIXED today
python3 testsys/perf/measure_ringing.py --case test.tpv37 --max-p2p 0.2
```

Variant binaries: copy `src/fortran` beside a `scripts/` symlink (the
makefile's `srcStamp` step runs `../../scripts/src_hash.py`), edit
`globalvar.f90:105` (`C_hg`) and `:191` (`kapa_hg`), `MACHINE=ubuntu make`.

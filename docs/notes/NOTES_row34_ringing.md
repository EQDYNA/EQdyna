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

## Follow-up 2026-09-29 (owner: "it seems to be related to wedge elements") -- controlled wedge vs hex A/B

**Verdict: the wedges are not the cause. Replacing tpv37's degenerate-wedge
mesh by the code's only hex alternative (planar dipping fault inserted into
hexes, `C_degen=0, insertFaultType=1`) at the same dip, dx, friction,
nucleation, stations and term makes the dp180 zigzag LARGER (rms 0.150 ->
0.197 MPa, p2p 0.673 -> 0.995) and moves the rupture-time field further from
the 50 m archive. Same picture at 30 deg. The mesh-scale ringing is a
property of the 500 m discretisation, not of the element type.**

### What had to change to make the hex meshing run at 15 deg

- tpv37 sets `dy = dx*cos(dip)` (483 m), which is exactly the fault's
  y-climb per down-dip node, so `case.setup` REFUSED the inserted geometry
  (`lib.FaultGeometryError: the fault surface climbs by 1 of a fault-normal
  cell ... >= 1.0 ... the elements inserted by insertFaultInterface tangle`)
  and the Fortran refused the same file with exit code 31. The comment in
  `tpv36_37_common.py` ("don't insert planar dipping fault because of the
  low dipping angle") is this refusal. With `dy = dx` (500 m) the ratio is
  cos(15) = 0.966 and the validator passes with its WARNING (|dy/dz| =
  3.73 > 0.2). Same at 30 deg (cos 30 = 0.866).
- At the case's `dt = 0.5*dz/vp` the inserted meshes went NaN at the fault's
  bottom edge (15 deg: node (500, 27046, -7247) m at step 215; 30 deg: node
  (500, 24249, -14000) m at step 208). The blend `ycoort = y*(ymax-peak)/ymax
  + peak` compresses the +y cells at the fault bottom to (1 - 27.0/37.0) =
  0.27 of dy; fault-normal thickness of fault-adjacent elements measured
  (Python port `meshgen.build_elements`, serial copy) min 48 m at 15 deg,
  94 m at 30 deg, against dz = 129 / 250 m. Runs completed at dt = 0.0025 s
  (15 deg, 2000 steps) and 0.005 s (30 deg, 1000 steps) -- 4.3x / 4.2x the
  wedge meshes' step count for the same 5 s.
- On-fault station FILES on the inserted mesh are named by VERTICAL depth
  (`library_output.f90:48-68`: with `C_degen==0` `fltxyz(2,4,1)` is 90 deg,
  so `dp180` (18 km down-dip) is written as `faultst000dp047`). Mapped by
  symlink for the comparison below; the station coordinates in
  `bStations.txt` are identical between the two meshings (diff clean).

### Element quality near the fault (Python port mesh, serial, fault-adjacent elements = any node in nsmp)

| mesh | dx, dy, dz (m) | fault-adjacent elems | of which wedges | scaled Jacobian at non-collapsed corners min / median | fault-normal thickness min / median / max (m) | edge-length ratio median |
|---|---|---|---|---|---|---|
| wedge 15 deg (tpv37 as committed) | 500, 483, 129 | 13900 | 6720 | 0.966 / 1.000 | 125 / 250 / 253 | 3.9 |
| hex-insert 15 deg, dy=dx | 500, 500, 129 | 7068 | 0 | 0.259 / 0.259 (= sin 15) | 48 / 132 / 494 | 1.7 (interior region max 16) |
| wedge 30 deg | 500, 433, 250 | 13900 | 6720 | 0.866 / 1.000 | 217 / 433 / 438 | 2.0 |
| hex-insert 30 deg, dy=dx | 500, 500, 250 | 7068 | 0 | 0.500 / 0.500 (= sin 30) | 94 / 254 / 882 | 1.7 |

The wedge mesh's own anisotropy is dz/dx = sin(dip) = 0.26 (129 m vertical
cells); the inserted hexes trade the collapsed corners for a 15 deg
parallelogram cross-section (scaled Jacobian sin(dip)) and an asymmetric
fault-normal thickness (48 m on +y, 494 m on -y at the fault bottom).
Neither meshing gives isotropic 500 m cells at the fault; an isotropic
control would need dx = dz along strike too (16x elements) and was not run.

### A/B at dip 15 (tpv37 physics), Fortran 4 ranks, 5 s, `v-shear-stress`

| station | metric | wedge (committed reference) | hex-insert (dy=dx, dt 0.0025) | 50 m archive |
|---|---|---|---|---|
| faultst000dp180 | ringing rms / p2p (MPa) | 0.150 / 0.673 | 0.197 / 0.995 | 0.0018 / 0.016 |
| | max\|run - archive\| (MPa) | 1.246 | 1.617 | -- |
| | t_arr (s) / peak slip rate (m/s) | 0.313 / 0.52 | 0.305 / 0.62 | 0.302 / 0.64 |
| faultst000dp120 | ringing rms / p2p | 0.061 / 0.247 | 0.155 / 0.691 | 0.024 / 0.123 |
| | t_arr | 3.127 | 2.635 | 3.087 |
| faultst040dp180 | ringing rms / p2p | 0.112 / 0.519 | 0.174 / 0.710 | 0.008 / 0.042 |
| | t_arr | 2.685 | 3.185 | 2.435 |
| faultst000dp240 | ringing rms / p2p | 0.161 / 0.900 | 0.495 / 2.030 | 0.009 / 0.041 |
| | t_arr | 3.354 | 2.920 | 3.045 |
| faultst080dp180 | t_arr | 4.918 | not reached by 5 s | 4.311 |
| whole fault (3477 nodes, matched to 1 m) | ruptured at 5 s / flips vs wedge | 1288 / -- | 1227 / 203 | -- |
| | \|d rupture time\| vs wedge mean / max (s) | -- | 0.330 / 1.017 | -- |
| | fault peak slip rate (m/s) | 5.01 | 6.61 (1.32x) | -- |

The inserted mesh is faster down-dip (dp120 2.64 vs archive 3.09 s) and
slower along strike (st040dp180 3.19 vs 2.44 s) -- a direction-dependent
distortion the wedge mesh does not have (its arrival errors vs the archive
are +0.04 / +0.25 / +0.31 s at dp120 / st040dp180 / dp240; the inserted
mesh's are -0.45 / +0.75 / -0.13 s).

### Mechanism control at dip 30 (NOT tpv37's physics; only the dip and hence dz change)

| station | wedge 30 rms / p2p (MPa), t_arr | hex-insert 30 rms / p2p, t_arr |
|---|---|---|
| faultst000dp180 | 0.114 / 0.556, 0.333 s | 0.183 / 0.837, 0.310 s |
| faultst000dp120 | 0.037 / 0.165, 3.271 s | 0.123 / 0.473, 2.740 s |
| faultst040dp180 | 0.048 / 0.255, 2.458 s | 0.136 / 0.574, 2.880 s |
| faultst000dp240 | 0.123 / 0.511, 3.354 s | 0.444 / 1.972, 2.965 s |
| whole fault | 1241 ruptured | 1226 ruptured, 157 flips, \|d rt\| mean 0.293 s, peak slip rate 1.34x |

Both dips, every station: the hex-only mesh rings 1.3-3x MORE than the wedge
mesh. Going from dip 15 to 30 (dz 129 -> 250 m, less anisotropic) lowers the
wedge mesh's dp180 rms by 24% (0.150 -> 0.114) -- consistent with the
resolution/anisotropy reading, and much smaller than the 3x from halving dx.

### What a fix would be

- Switching tpv36/37 to `insertFaultType=1` is not a fix: it needs `dy=dx`
  and a 4x smaller dt to run, rings more, and distorts the rupture-time
  field more than the wedges do. It would also change both cases' meshes and
  every committed reference (rule 7, owner's call) for a worse result.
- A wedge-specific hourglass treatment (Fortran + port, rule 23) has no
  target: the ringing survives, larger, on a mesh with zero wedges.
- What does move it is dx (3x per halving; 50 m is flat). Fault-normal
  anisotropy (dz = dx sin dip) is the remaining untested suspect; the
  controlled test is a wedge mesh with dx = dz along strike (~16x cost at
  dip 15) and was not run.

### Reproduce (scratch, this branch)

Case copies were made with `create.newcase <name> test.tpv37` and `sed` on
the copied `tpv36_37_common.py`: `par.C_degen = 0`, `par.insertFaultType =
1`, `par.fymin, par.fymax = 0.0, 0.0`, `par.dy = par.dx`, `par.dt = 0.0025`
(15 deg) / `0.005` (30 deg); `par.dip = 30` for the 30 deg pair; `par.nx =
par.nz = 1` copies for the Python-port mesh-quality read. Ringing via
`testsys/perf/measure_ringing.py --stations ...` on a directory of symlinks
mapping `faultst000dp047 -> faultst000dp180` etc.; rupture fields matched
node-by-node (rounded to 1 m) between canonicalised frt sets.

## C_hg=3: KF78 + Flanagan-Belytschko viscous hourglass control (owner: "make FB chg=3", 2026-09-29)

Implemented on this branch, entirely inside `src/fortran/calcHourglassResist.f90`
(no change to `assembleGlobalMass.f90` / `calcSSPhi4Hrgls`; `globalvar.f90:105`
comment only). C_hg=1 and 2 arithmetically untouched: C_hg=3 with kapa_hg=0
reproduces the committed tpv37 reference at max|diff| = 0.0 over all 22 frt
columns and 0.0000 MPa at dp180; unit + regression tiers green.

- Vectors: `calcFBHourglassVectors` (end of file) computes, per element,
  gamma_aI = Gamma_aI - sum_i (sum_J Gamma_aJ x_iJ) dN_I/dx_i (FB 1981, IJNME
  17:679-706, eq. 3.44) from `meshCoor` and the stored one-point gradients
  `eleshp(1:3,:,nel)`, then RMS-normalises each gamma_a so sum_I gamma_aI^2 = 8
  = sum_I Gamma_aI^2 (the normalisation calcSSPhi4Hrgls applies to KF78's phi;
  it gives gamma the raw Gamma's norm so the C_hg=2 coefficient applies
  unchanged). On a rectangular brick gamma == Gamma.
- Force (C_hg==3 block): q_ia = sum_I v_iI gamma_aI, f_iI = -c sum_a q_ia
  gamma_aI, c = kapa_hg*rho*Vp*V^(2/3)/4 with V = eledet*w (LS-DYNA hourglass
  type 2 / Goudreau-Hallquist 1982 coefficient), applied to the velocity `vl`,
  assembled with the same sign convention as the C_hg=2 branch.
- Exposure: C_hg and kapa_hg remain compile-time defaults in `globalvar.f90`;
  runs below used scratch copies with the two lines edited.

### tpv37, 500 m, Fortran 4 ranks, 5 s, down-dip shear (rms / p2p MPa; ref = committed C_hg=1 reference)

| kapa_hg | dp180 | dp120 | dp240 | st040dp180 | max\|dp180 - 50 m archive\| | ruptured (ref 1288) / flips | mean / max \|d rt\| s | peak slip rate ratio |
|---|---|---|---|---|---|---|---|---|
| 0 (= C_hg=1) | 0.150 / 0.673 | 0.061 / 0.247 | 0.161 / 0.900 | 0.112 / 0.519 | 1.246 | 1288 / 0 | 0 / 0 | 1.000 |
| 0.03 | 0.082 / 0.531 | 0.041 / 0.155 | 0.130 / 0.649 | 0.093 / 0.525 | 1.254 | 1270 / 18 | 0.014 / 0.205 | 0.897 |
| 0.05 | 0.065 / 0.449 | 0.038 / 0.147 | 0.122 / 0.574 | 0.090 / 0.496 | 1.258 | 1258 / 30 | 0.022 / 0.205 | 0.852 |
| 0.07 | 0.053 / 0.377 | 0.035 / 0.142 | 0.115 / 0.521 | 0.086 / 0.470 | 1.259 | 1256 / 32 | 0.029 / 0.216 | 0.814 |
| 0.10 | unstable (46 Hz, 55 MPa p2p at dp120, peak slip rate 5.7e4 m/s) | | | | | 3477 / 2189 | | 1.1e4 |
| 0.15 | blows up (1e10 MPa) | | | | | | | |
| 50 m archive | 0.0018 / 0.016 | 0.024 / 0.123 | 0.009 / 0.041 | 0.008 / 0.042 | 0 | | | |

The geometry-corrected gamma changes essentially nothing relative to the
earlier scratch raw-Gamma "KF78+viscous" run: at kapa 0.03, rms 0.0822 vs
0.0828, 18 flips both, peak-slip-rate ratio 0.897 vs 0.886; at kapa 0.1 the
same 46 Hz instability at the same amplitude. Stable range on this mesh at
its CFL 0.5: kapa_hg <= 0.07 (0.1 unstable). Corrected vectors do not move
the stability limit.

### Clean-case control: test.tpv8 (hex, vertical) at C_hg=3, kapa_hg=0.05, vs its committed reference

1891 nodes: ruptured 830 -> 767, **63 flips**, mean |d rupture time| 0.103 s
(max 0.35 s), fault peak slip rate ratio 0.987; frt max|diff| 2.7e7 against
CASE_BOUND 1e-8 (i.e. the gate would fail, as any physics change must).
faultst000dp075 h-shear detrended rms 0.723 -> 0.217 MPa over [0.59, 2.59] s:
the term damps post-front content on a clean hex case too, and shifts its
rupture arrivals by 0.1 s on average.

### Verdict

C_hg=3 works as designed and is stable to kapa_hg 0.07, but it does NOT kill
the zigzag and it DOES move the rupture: at the largest stable coefficient
the dp180 ringing rms is still 29x the 50 m archive's, the archive error at
dp180 is unchanged (1.25 -> 1.26 MPa: it is front-shape, not ringing), and
the price is 32 rupture-flag flips, -19% fault peak slip rate on tpv37 and
63 flips / 0.10 s mean arrival shift on tpv8. In the time series
(archive / default / kapa 0.05 at four stations, plot not kept) the damped
curve is the default one at ~60% amplitude, same frequency, not the flat
archive curve.

### What a PR would need, if the owner still wants the option landed

- Fortran: this file as is (`calcHourglassResist.f90`), plus reading
  `C_hg`/`kapa_hg` from input (`readInputFiles.f90`, `bGlobal.txt` line
  order), `scripts/case.setup` writer and `scripts/defaultParameters.py`
  defaults (C_hg=1, kapa_hg=0.1 today are compile-time).
- Python port (rule 23): `src/python/eqdyna/assembleGlobalKU.py::
  calcHourglassResist` is KF78-only; needs the gamma computation (vectorised
  over elements from `meshCoor[conn]` and `eleshp`) and the viscous scatter
  as a second contribution in the fused loop, with the correspondence table
  entry updated; parity gate = C_hg=3 fortran vs python-jax on tpv37.
- A regression guard that C_hg=3 with kapa_hg=0 is bit-identical to C_hg=1
  (the check done by hand above), and one that C_hg=1/2 outputs are unchanged.
- Adoption is the owner's call; every adopting case changes physics and needs
  a new committed reference (rule 7). On this evidence no gated case should
  adopt it.

## Closed 2026-09-30 (owner: "just use 500 m as regression as it is")

Row 34 closed: tpv36/37 stay gated at 500 m with the ringing accepted;
refinement is what removes it (table above). The local research branch
`research/ringing-36-37` (tip `515db32`, never pushed) was deleted after this
record was written. Its only code not on master is C_hg=3, kept here verbatim
so it can be re-applied with `git apply` from a checkout of `5647889`:

```diff
diff --git a/src/fortran/calcHourglassResist.f90 b/src/fortran/calcHourglassResist.f90
index 1a2a62f..409863f 100644
--- a/src/fortran/calcHourglassResist.f90
+++ b/src/fortran/calcHourglassResist.f90
@@ -7,7 +7,7 @@ subroutine calcHourglassResist
     include 'mpif.h'
 
     integer (kind = 4) :: nel, i , j, k, itmp, itag, fi(4,8)
-    real (kind = dp) :: phid(ned), dl(ned,nen), vl(ned,nen), fhr(ned,nen), f(24), det, coef, q(3,4)
+    real (kind = dp) :: phid(ned), dl(ned,nen), vl(ned,nen), fhr(ned,nen), f(24), det, coef, q(3,4), fvis, gam(nen,4)
     
     ! LOCAL stage timer. The shared global startTimeStamp was written by
     ! eight sites across six files, each pairing it with a different
@@ -21,7 +21,11 @@ subroutine calcHourglassResist
                 dl(j,i) = dispArr(j,nodeElemIdRelation(i,nel)) + rdampk* vl(j,i)
             enddo
         enddo
-        if (C_hg == 1) then
+        ! C_hg == 3 (Flanagan-Belytschko viscous, 2026-09-29, row 34) runs the
+        ! KF78 stiffness branch below UNCHANGED and then adds the viscous
+        ! block at the end of this element loop. C_hg == 1 and 2 are
+        ! arithmetically untouched by this addition.
+        if (C_hg == 1 .or. C_hg == 3) then
             do itmp = 1, 4
                 fhr = 0.0d0
                 !... calculate sum(phi*dl)
@@ -60,6 +64,15 @@ subroutine calcHourglassResist
 
         elseif (C_hg==2) then
             !viscous hourglass control
+            ! MEASURED 2026-09-29 (docs/notes/NOTES_row34_ringing.md): NOT
+            ! recommended for test.tpv36/test.tpv37 or any distorted or
+            ! degenerate-wedge mesh. This branch REPLACES the KF78 stiffness
+            ! control and uses the uncorrected +-1 Gamma vectors; on tpv37 at
+            ! 500 m it collapses the rupture at every kapa_hg tried,
+            ! including kapa_hg = 0 (1288 -> 614/541/520/436 ruptured nodes
+            ! at kapa 0/0.01/0.03/0.1), and rings MORE than C_hg=1 at
+            ! kapa_hg <= 0.03. Use C_hg = 3 for a viscous term that keeps
+            ! KF78 and uses the geometry-corrected vectors.
             coef = 0.25d0*kapa_hg*mat(nel,3)*mat(nel,1)*(eledet(nel)*w)**(2.0d0/3.0d0)
             fi(1,1)=1;fi(1,2)=1;fi(1,3)=-1;fi(1,4)=-1;fi(1,5)=-1;fi(1,6)=-1;fi(1,7)=1;fi(1,8)=1
             fi(2,1)=1;fi(2,2)=-1;fi(2,3)=-1;fi(2,4)=1;fi(2,5)=-1;fi(2,6)=1;fi(2,7)=1;fi(2,8)=-1
@@ -95,6 +108,104 @@ subroutine calcHourglassResist
                 enddo
             enddo
         endif
+
+        if (C_hg == 3) then
+            ! Flanagan & Belytschko (1981, IJNME 17:679-706) viscous hourglass
+            ! control, added on top of KF78 (above). The geometry-corrected
+            ! hourglass shape vectors gamma are computed HERE, per element,
+            ! by calcFBHourglassVectors (end of this file) from the element's
+            ! node coordinates and its stored one-point shape-function
+            ! gradients eleshp -- nothing is added to assembleGlobalMass.f90
+            ! or calcSSPhi4Hrgls. Coefficient (LS-DYNA theory manual,
+            ! hourglass type 2 "Flanagan-Belytschko viscous form"; Goudreau
+            ! & Hallquist 1982), identical to the C_hg == 2 branch's:
+            !   q_ia  = sum_I v_iI gamma_aI          (hourglass velocity rates)
+            !   f_iI  = - c * sum_a q_ia gamma_aI
+            !   c     = kapa_hg * rho * Vp * V^(2/3) / 4,  V = eledet*w
+            ! Acts on the VELOCITY vl (not on dl = disp + rdampk*vel, which
+            ! is the KF78 stiffness argument). For a rectangular brick gamma
+            ! == Gamma and this block equals the C_hg == 2 force; on a
+            ! distorted or degenerate (wedge) element gamma is orthogonal to
+            ! the linear velocity field where Gamma is not, so it damps only
+            ! the zero-energy modes.
+            call calcFBHourglassVectors(nel, gam)
+            coef = 0.25d0*kapa_hg*mat(nel,3)*mat(nel,1)*(eledet(nel)*w)**(2.0d0/3.0d0)
+            q = 0.0d0
+            do itmp = 1, 4
+                do i = 1, ned
+                    do j = 1, nen
+                        q(i,itmp) = q(i,itmp) + vl(i,j)*gam(j,itmp)
+                    enddo
+                enddo
+            enddo
+            do i = 1, nen
+                do j = 1, ned
+                    fvis = 0.0d0
+                    do itmp = 1, 4
+                        fvis = fvis - coef*q(j,itmp)*gam(i,itmp)
+                    enddo
+                    if (numOfDofPerNodeArr(nodeElemIdRelation(i,nel)) == 3) then
+                        itag = eqNumStartIndexLoc(nodeElemIdRelation(i,nel))+j
+                    elseif (numOfDofPerNodeArr(nodeElemIdRelation(i,nel)) == 12) then
+                        itag = eqNumStartIndexLoc(nodeElemIdRelation(i,nel))+j+9
+                    endif
+                    k = eqNumIndexArr(itag)
+                    if(k > 0) then
+                        nodalForceArr(k) = nodalForceArr(k) + fvis
+                    endif
+                enddo
+            enddo
+        endif
     enddo
     compTimeInSeconds(5) = compTimeInSeconds(5) + MPI_WTIME() - tStageStart
 end subroutine calcHourglassResist
+
+
+subroutine calcFBHourglassVectors(nel, gam)
+    ! Flanagan & Belytschko (1981, IJNME 17:679-706) hourglass shape vectors
+    ! for the one-point 8-node hexahedron:
+    !   gamma_aI = Gamma_aI - sum_i ( sum_J Gamma_aJ x_iJ ) dN_I/dx_i     (FB81 eq. 3.44)
+    ! with Gamma the four +-1 hourglass base vectors (FB81 Table 2, the same
+    ! ordering as calcSSPhi4Hrgls's `ha` and the C_hg==2 branch's `fi`), x_iJ
+    ! the element node coordinates and dN_I/dx_i the one-point shape-function
+    ! gradients this element already stores in eleshp (assembleGlobalMass.f90,
+    ! assembleElementMassDetShg). gamma is orthogonal to every linear
+    ! velocity field on the actual (distorted or degenerate) element, which
+    ! the raw Gamma is only on a parallelepiped.
+    !
+    ! NORMALISATION (stated because the viscous coefficient depends on it):
+    ! each gamma_a is RMS-normalised, gamma_a <- gamma_a / sqrt(sum_I
+    ! gamma_aI^2 / 8), so sum_I gamma_aI^2 = 8 = sum_I Gamma_aI^2. This is the
+    ! same normalisation calcSSPhi4Hrgls applies to KF78's phi, and it makes
+    ! gamma carry exactly the norm of the raw Gamma, so the Goudreau-Hallquist
+    ! coefficient c = kapa_hg*rho*Vp*V^(2/3)/4 of the C_hg==2 branch applies
+    ! unchanged. On a rectangular brick gamma == Gamma identically.
+    use globalvar
+    implicit none
+    integer (kind = 4), intent(in) :: nel
+    real (kind = dp), intent(out) :: gam(nen,4)
+    integer (kind = 4) :: a, i, j, k
+    integer (kind = 4), dimension(8,4) :: ha = reshape((/ &
+            1,1,-1,-1,-1,-1,1,1, 1,-1,-1,1,-1,1,1,-1, &
+            1,-1,1,-1,1,-1,1,-1, -1,1,-1,1,1,-1,1,-1/), &
+            (/8,4/))
+    real (kind = dp) :: gx(3), rms
+
+    do a = 1, 4
+        gx = 0.0d0
+        do j = 1, nen
+            do i = 1, 3
+                gx(i) = gx(i) + ha(j,a)*meshCoor(i,nodeElemIdRelation(j,nel))
+            enddo
+        enddo
+        rms = 0.0d0
+        do k = 1, nen
+            gam(k,a) = ha(k,a) - (gx(1)*eleshp(1,k,nel) + gx(2)*eleshp(2,k,nel) + gx(3)*eleshp(3,k,nel))
+            rms = rms + gam(k,a)**2
+        enddo
+        rms = sqrt(rms/8.0d0)
+        do k = 1, nen
+            gam(k,a) = gam(k,a)/rms
+        enddo
+    enddo
+end subroutine calcFBHourglassVectors
diff --git a/src/fortran/globalvar.f90 b/src/fortran/globalvar.f90
index 5962330..8860cec 100644
--- a/src/fortran/globalvar.f90
+++ b/src/fortran/globalvar.f90
@@ -102,7 +102,7 @@ MODULE globalvar
     integer (kind = 4) :: C_elastic              ! 1 = elastic version; 0 = plastic version
     integer (kind = 4) :: C_nuclea                ! 1 = allow artificial nucleation; 0 = disabled
     integer (kind = 4) :: C_Q  = 0                ! only with C_elastic==1: 1 = allow Q attenuation; 0 = do not
-    integer (kind = 4) :: C_hg = 1                ! hourglass control: 1 = KF78, 2 = viscous HG
+    integer (kind = 4) :: C_hg = 1                ! hourglass control: 1 = KF78, 2 = viscous HG (raw Gamma, replaces KF78), 3 = KF78 + Flanagan-Belytschko viscous (row 34)
     integer (kind = 4) :: C_dc = 0                ! double-couple source: 1 = yes, 0 = no
     integer (kind = 4) :: C_degen                  ! degenerate-element flag: 0 = brick, 1 = wedge, 2 = tetra
     integer (kind = 4) :: output_plastic          ! 1 = write plastic-strain output
```

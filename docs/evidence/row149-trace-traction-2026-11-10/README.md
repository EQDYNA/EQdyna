# Row 149: TPV12/13 free-surface-trace normal stress, research attempt of 2026-11-10

Question: why does the normal stress at `faultst000dp000` (the node on both the
60-degree dipping fault and the free surface) diverge with grid refinement, and
is there a physically derived trace-node treatment that stops it?

## Verdict: no fix. A real inconsistency was found and repaired in a scratch build; the outcome curve did not move (rule 4a)

Nothing in `src/` changed. No reference was touched. No 100 m run was made,
because the condition for one (the coarse curve converging) was not met.

## Runs

EQdyna fortran, 4 MPI ranks, sources of master `6559464`, each run to 4.5 s,
one run per cell (run-to-run spread not measured). `dp000_nstress.csv` holds
every run; the values below are n-stress at `faultst000dp000` at 4 s, in MPa.

| dx (m) | TPV13 master | TPV13 spike | TPV12 master |
|---|---|---|---|
| 500 | -11.7612 | -13.6868 | -9.2340 |
| 250 | -2.5297 | not run | not run |
| 200 | +1.9555 | -3.4215 | -3.9865 |
| 150 | +5.4215 | -0.2087 | -2.2977 |

FaultMod at 100 m, same station, 4 s: TPV13 +1.5193, TPV12 -3.0050 (8 s:
+2.9261, -2.3820). It is the only other-code submission in the shared dataset
and exists at one resolution, so its own behaviour under refinement is not
measurable here. Its file header gives the station as "on-fault, 0 m along
strike, 0 m depth" and does not say whether the value is nodal or a face
average.

At depth the two codes agree: `faultst000dp015` at 4 s is -26.41 (EQdyna 150 m)
against -27.4491 (FaultMod 100 m).

## What was checked and found sound

- The split-node traction is the constraint force over the node's tributary
  area (`src/fortran/faulting.f90`), and the trace node gets the half-row area
  (`src/fortran/meshgen.f90`, `arn`): 11250 m^2 at the trace against 22500 m^2
  one row down at 150 m.
- The instrumented build (`diag_trace_instrumentation.patch`) reproduces the
  master build's station value to every printed digit at 500 m and 200 m.

## Where the tension comes from

The instrumentation splits the trace-node normal traction into the element
stress force (KU) and the hourglass-resistance force (HG). TPV13, 4 s, MPa:

| dx (m) | KU | HG | total |
|---|---|---|---|
| 250 | -9.190 | +6.660 | -2.530 |
| 200 | -6.729 | +8.684 | +1.955 |
| 150 | -4.636 | +10.057 | +5.422 |

The elastic case (TPV12) carries less than half of that hourglass part at the
same dx (+3.646 at 200 m, +4.010 at 150 m).

## The hypothesis, the spike and its falsification

Hypothesis: `calcElemKU.f90` returns the one-point deviatoric stress to the
Drucker-Prager surface, while `calcHourglassResist.f90` keeps charging the
elastic stiffness on total displacement, so the non-constant part of the
element stress never yields and the trace tension grows with refinement.

Spike (`spike_plastic_hourglass.patch`, plastic cases only): the hourglass
generalized displacement is integrated incrementally and relaxed each step by
the element's own return factor; with no yielding it reduces to the existing
form.

Kill criterion, written down before the 200 m and 150 m results: the 200 m to
150 m rise of dp000 must shrink from the master build's +3.47 MPa to under
about 1 MPa.

Result: the spike shifts dp000 by -5.4 to -5.6 MPa at 200 m and 150 m, but the
200 m to 150 m rise is +3.21 MPa. The defect is real; it is not the cause of
the divergence. The elastic control agrees: TPV12 has no plasticity and its
trace value still rises with refinement (-9.23, -3.99, -2.30).

## What does converge

Mean normal traction over a fixed down-dip length from the trace (piecewise
linear through the nodal values, `ana.py`), TPV13 master, 4 s, MPa:

| averaging length | 250 m grid | 200 m grid | 150 m grid |
|---|---|---|---|
| top 250 m | -9.931 | -8.737 | -8.075 |
| top 500 m | -15.375 | -15.071 | -14.869 |
| top 1000 m | -18.801 | -18.690 | -18.514 |

The jump between the top two fault nodes stays near 15-18 MPa at every dx in
TPV13: the steep gradient sits inside the top element at any resolution, and
the trace node samples an ever-smaller piece of it. The resultant over 500 m
or more is stable to about 3% across these grids; the nodal value is not a
convergent quantity.

## Side finding, not traced to a line

TPV12's trace node starts with zero normal traction (0.0000 at every dx).
TPV13's starts at -1.8410 / -0.7379 / -0.5536 at 500 / 200 / 150 m, which
scales to -0.369 at 100 m, FaultMod's initial value (-0.3695) for both cases.

## Not run

- No python-jax port, regression guard, unit/regression tier or e2e cells:
  there is no treatment to gate.
- No 100 m run, no 8 s values, no frt comparison of the spike against master.
- The 100 m point quoted on the board (+8.43 MPa) is the October run and was
  not re-measured.

## Reproduce

Apply `diag_trace_instrumentation.patch` to `src/fortran/` in a scratch copy,
build, run `test.tpv13` at the dx wanted with `par.term = 4.5`, then
`python3 ana.py <run_dir> 4.0`. Each patch is a diff against master `6559464`
on its own. Neither patch is a proposal to merge.

# Session log — 2026-09-22, the release session (v5.15.0, v5.13.0 backfill, v5.16.0)

Conductor: wei-lin. Budget granted 20:08, 12 hours, standing rolling renewal.
Session opened 21:43 at `d7fe2c2` + two dirty files; publish authority granted
by the owner for v5.15.0 and v5.16.0 specifically.

Grant, stated back before the first tag: minor and patch tags on `master` of
this repo, plus `gh release create` for v5.15.0, v5.16.0 and the retroactive
v5.13.0 object. No major bump, no force-update or re-point of an existing tag
(rule 8), nothing outward-facing beyond this repo.

---

## 0. The two dirty files — landed first, unedited

`PROJECT_RULES.md` and `pathway_forward.md` were dirty in the main checkout,
authored by zofia-kaminska (rule 15b; board item 69) and left uncommitted.
Uncommitted rule text is the exact hazard rule 20b names, and here it was the
rules file itself.

Committed unedited at **`87af56e`**, pushed. Count verified against the file's
own self-check before committing: `grep -c '^## ' PROJECT_RULES.md` = 31,
`grep -c '^## [0-9]*\. '` = 21 — both match what the amended text claims.

**Rule 15b's workflow change is NOT landed.** The rule text is in tree; the
`test.yml` edit it describes is not. It stays text-only until its two guards
exist and until the owner sets the staleness bound for the ancestor CI run
(zofia proposed 7 days, labelled it a guess, and said treat it as 0 until he
rules). An operator-chosen reduced CI run without those guards is the
silently-skipped step rule 2 forbids. The rule's own routing clause agrees:
the workflow edit goes to whoever conducts the release, and this conductor
declines it pending the gate.

---

## 1. DEVIATION FROM THE HANDOFF — the rough-fault fix was already on master

Recorded loudly, per the loop rule that code beats plan when the code is
unambiguous. The handoff described v5.16.0 as "the rough-fault normal fix
alone (`origin/wei/rough-fault-normal-2026-09-22`)", still to be landed, and
v5.15.0 as "89 commits from master HEAD, NO reference moves". Both halves of
that picture were stale. What the tree actually says:

| fact | evidence |
|---|---|
| the v5.15.0 release commit already EXISTS on master | `54e0697` "release: v5.15.0 -- GPU jax backend, the perf ledger, and a pre-tag CI gate"; README + VERSION + pastReleaseNotes + the Fortran banner, all four, in one commit |
| it was already CI-green and merely untagged | session log `12c5be0` "v5.15.0 CI-green and untagged"; re-verified by me, below |
| the rough-fault normal fix is ALSO already on master | `44afd4a` "Merge: rough-fault normals recomputed at the mesh dx", `git merge-base --is-ancestor 44afd4a HEAD` = yes |
| it carries the two reference regenerations | `bbe6f1a` (tpv29), `642f119` (tpv30), plus `51d37a0` for the shipped surfaces |
| `54e0697` predates all of it | `git merge-base --is-ancestor 54e0697 702b56b` = yes |
| v5.14.0..`54e0697` moves NO reference and NO case input | `git diff --name-only v5.14.0..54e0697 \| grep -E 'test.reference.results\|case_input'` → empty |
| the count is 35, not 89 | `git rev-list --count v5.14.0..54e0697` = 35; `..HEAD` = 90 |

**The owner's intent survives intact and is unchanged; only the mechanism
changed.** He asked for two tags so that a reference regeneration is not
buried in a release of unrelated commits — "if someone later asks when tpv29's
reference moved and why, the answer must be one tag with one reason." That is
achieved by tagging v5.15.0 at `54e0697` (the pre-rough-fault release commit,
no reference moves) rather than at master HEAD, and cutting v5.16.0 from HEAD.
A tag points where it is told; nothing had to be unwound, rebased or
re-landed, and no gate was skipped to save time.

What I did NOT do, and will not: fold the two into one release. The split is
the instruction and it is also right.

---

## 2. v5.15.0 — tagged at `54e0697`, published 21:50

Rule 15a run in full, BEFORE `git tag`, against the exact SHA:

    python3 testsys/regression/check_pretag_ci.py --pre-tag 54e0697
    PASS  a completed, successful Automatic Testing of EQdyna run exists for
          54e0697865ac798d8c6d5dc90c335af92e9b071b
    exit 0

No `--ack-paths-ignored-parent` was needed or used: `54e0697` touches
`VERSION` and `src/fortran/eqdyna3d.f90`, neither in `test.yml`'s
`paths-ignore`, so it triggered and completed a CI run of its own.

Tag + push + `gh release create --verify-tag --latest` run as ONE
uninterrupted action, no pause to watch CI between them — the v5.8.2
`check_network_side` window (rule 15 step 7) closes only if they are chained.
Release body taken verbatim from the README News block at `54e0697`, minus the
pointer line.

<https://github.com/EQDYNA/EQdyna/releases/tag/v5.15.0>

Gate evidence for the tree itself is the release commit's own sweep, recorded
in its notes: 31 of 31 cells SUCCESS, exit 0, 3839.6 s, on the exact tree the
tag contains, with no source file changed after it. I did not re-run that
sweep: `54e0697` is an immutable historical commit whose gate run and CI run
both exist and both completed, and re-running a sweep against a tree that
cannot have changed buys nothing. The sweep I *am* running (below) is for
v5.16.0, whose tree has genuinely moved.

---

## 3. v5.13.0 — the missing Release object, backfilled (board item 54)

`v5.13.0` was tagged 2026-09-19 and never got a GitHub Release, so the
Releases page skipped from v5.12.0 to v5.13.1. Backfilled from the retroactive
notes already written into `pastReleaseNotes.md`:

    gh release create v5.13.0 --notes-file <block> --latest=false --verify-tag \
      --title "v5.13.0 (retroactive notes)"

<https://github.com/EQDYNA/EQdyna/releases/tag/v5.13.0>

`--latest=false` explicitly, so backfilling an old tag cannot demote v5.15.0 on
the Releases page. The tag itself was NOT re-pointed and its tree stays
permanently self-inconsistent (it declares 5.12.0 at runtime) — that is rule 8
and it is recorded in `pastReleaseNotes.md` rather than hidden.

---

## 4. v5.16.0 — in flight

Contents: the rough-fault normal fix (`44afd4a`), carrying the tpv29 and tpv30
reference regenerations, plus the docs/figures/board commits after it.
`git diff --name-only 44afd4a..HEAD` touches no solver source, no reference and
no `testsys/` file — only `docs/`, `scripts/figures/`, the board and the rules.
So master HEAD's gated surface is `44afd4a`'s.

Headline, from the handoff and to be restated in the notes only after my own
sweep confirms the cells: at dx=100 both cases reproduce the owner's 2015 SCEC
submissions — tpv29 median |dt| 0.0073 s, tpv30 0.0107 s, both inside the
archive's own 100 m-vs-50 m sensitivity of 0.0160 s; rupture-time seeds 28 → 3.
Cause was derivative columns subsampled rather than recomputed at the mesh
spacing.

Gate: a fresh full `run.py all` on `87af56e`, launched detached 21:48
(rule 20), polled by its log artifact rather than by its process. A preserved
gate run for `44afd4a` exists at `docs/evidence/gate-44afd4a/`, and I am not
landing a release on it — that is a report, and rule 15 step 1 asks for a green
gate on the tree being tagged, which is mine to run.

Sequence from here, no step skippable: sweep green → zofia writes the
Tasks-done row and board rows → release commit (VERSION 5.16.0, README News,
pastReleaseNotes move, Fortran banner) pushed ALONE → CI green on that exact
SHA → `check_pretag_ci.py` → tag + push + `gh release create` as one action.

If the budget runs out between v5.15.0 and v5.16.0, the instruction is to stop
after v5.15.0. It has not.

---

## 5. Dispatched during the sweep window

- **lars-eriksson**, read-only, no worktree (he holds no Edit/Write tools, and
  a worktree would have put the count at 5 against the owner's ceiling of 4):
  does `test.drv.a6`'s FRACTAL roughness path (`insertFaultType=2`) carry the
  same class of defect as the one `702b56b` fixed — derivative columns
  inconsistent with the mesh dx at which the solver consumes them? This is the
  owner's named top unaddressed risk. Briefed with the shared reader
  (`readInputFiles.f90:228-360`), both consumers, the generator candidates, the
  new guard, and board item 24(a)'s unexplained "p10-p90 spread 0.93-1.08,
  fractal roughness perturbing the local fault normal" — with the explicit
  warning that the earlier attribution of that symptom to the derivative bug
  was RETRACTED, so the link must be established or refuted from the code, not
  assumed either way. Told plainly that a real finding moves drv.a6's flip
  budget of 450 and its reference, which is a rule 7 reviewed change and not
  something to slip in.
  Constrained to no simulation runs: a sweep is in flight and item 42 forbids
  stacking two heavy jobs.

Not dispatched, deliberately: the jax-CPU perf A/B on
`mira/jaxmpi-step-2026-09-22`, because the owner's standing condition is an
idle box and the sweep owns it. GPU is off the queue entirely by his
instruction.

---

## 6. Owner-gated, decided by nobody here

Items 63, 56, 57, 17, 19(b), 32, and rule 15b's staleness bound. Untouched.

---

## 7. THE INCIDENT — a 40-minute gate destroyed by a concurrent session

The first v5.16.0 gate sweep was launched at 21:48 **in the main checkout**.
That was my error and rule 21a now exists because of it.

`testsys/e2e/run_e2e.py:587-594` rotates `$REPO_ROOT/test` to `test.prev` at
startup, unconditionally, with no lock of any kind. At **22:12:20** a second
Claude session ran its own e2e in the same checkout — its artifact,
`docs/perf_snapshots/e2e_cells_2026-09-22_221220_2462614.json`, landed on
master in `88a4012` — and rotated my in-flight tree out from under the sweep.

| | |
|---|---|
| cells already passed | 25 |
| cell killed by the rotation | `test.tpv1053d x python-numpy`, at **1500.1 s** |
| error | `FileNotFoundError: .../test/test.tpv1053d.python-numpy/frt.txt0` |
| raised at | `src/python/eqdyna/library_output.py:120`, from `eqdyna3d.py:522` — it writes by ABSOLUTE path, so the rename broke it at the end of the run |
| cells doomed but still running | 4 (`tpv29`, `tpv36`, `tpv37`, `drv.a6`, all python-numpy) |
| disposition | killed by PID after a targeted `ps`; never `pkill -f` |

**Why this is worse than a lost half-hour.** The failure surfaces as
`FAIL test.tpv1053d x python-numpy` in the sweep summary — visually identical
to a genuine parity regression. Had I taken it at face value I would have
reverted a good change or re-gated a landing that was already fine. The proof
it was infrastructure is the same cell, same SHA, in isolation twenty minutes
later: **SUCCESS in 907.0 s, max|diff| 4.569608e-06 against bound 1.0e-04.**

Both logs are kept under `docs/evidence/gate-v5.16.0/` (rule 8) — the passing
one and the collided one, the latter named `_INCIDENT` so it cannot be mistaken
for a gate.

Diagnosis discipline worth recording: my first two hypotheses were both wrong
and both were killed by evidence rather than by argument. "A second sweep is
running" was refuted by `ps` showing exactly one `run_e2e.py`. "The harness
cleans up cell directories" was refuted by reading the harness — the only
`rmtree`/`move` pair in it is the startup rotation. The actual cause only
appeared when I looked at *who else* had written to the repo, and the other
session's own committed snapshot timestamped it to the second.

Follow-ups, neither done here: pathway item 70 (a lock on `$REPO_ROOT/test`,
routed to `iris-vermeulen`) and rule 21a, which is hortatory until item 70
lands and says so about itself.

---

## 8. The gate that counts

Re-run in an ISOLATED worktree, `.claude/worktrees/wei-gate-v5160`, detached at
`88a4012`:

    ./install-eqdyna.sh -m ubuntu        exit 0
    python3 testsys/run.py all
      SUCCESS unit / SUCCESS regression / SUCCESS e2e
      ran: 31 of 40 cells in the 10 case x 4 backend table (31 passed, 0 failed)
    SWEEP EXIT: 0        22:25:50 -> 22:55:51  (1801 s)

31/31 is the full count rule 15 step 1 requires
(`len(matrix.CASE_BOUND) * 3 + len(matrix.PY_MPI_RANKS)`). The 31 cell
wall-clock rows and the snapshot were carried out of the worktree into the
release commit (rule 20b); ledger 224 -> 255 lines.

`88a4012..HEAD` touches only the board, the rules, `docs/` and the release
files — no solver source, no reference, no `testsys/` file. The gated code
surface IS the tagged code surface.

---

## 9. Two subagent findings, both re-verified by me at file:line

**The owner's named top unaddressed risk is CLOSED, refuted from the code.**
`test.drv.a6`'s fractal roughness path does not carry the `702b56b` defect
class and structurally cannot: `scripts/generateFaultInterface:94-95`
differentiates the surface at the same spacing it was generated on, so there is
no finer source to subsample from. Measured residual against a recomputation at
the mesh dx: **1.97e-16** (vs the guard's `STORED_DERIV_ATOL = 1.0e-4` and the
3.9e-3 the TPV29 bug actually moved those columns by). The covering check is
`scripts/lib.py:816-817` via `case.setup:379-380`, which runs for every
`insertFaultType > 0`, generated files included. Item 24(a)'s 0.93-1.08 spread
is genuine mesh-dependent fractal roughness, not a normal inconsistency — the
link to the TPV29 defect is refuted, and the earlier retraction was correct.
No budget and no reference moved. I re-read all four load-bearing lines myself
before recording it.

**The setup-rewrite premise is wrong, and this is a rule-2-of-the-loop
deviation I am recording rather than escalating.** The owner's direction —
"each rank builds only its own subdomain, following Fortran; `mesh4num.f90`
counts locally, `meshgen.f90` builds only that slab" — describes what the
Fortran **already does**. Verified directly, not relayed:

- `src/fortran/meshgen.f90:28-31` slices per-rank via `calcXyzMPIId`; the build
  loop at `:65-67` is `do ix = 1, nx` with `nx` local, sized at `:543-551` from
  `MPIXyzId`.
- `src/fortran/countMeshEntities.f90` — the routine the owner calls
  `mesh4num.f90`; it was **renamed** — does the identical local slicing at
  `:13-16` and assigns from purely local counters at `:72-75`.
- `src/fortran/library.f90:19-20` states it outright: "`totalNumOfElements` is
  THIS rank's count."
- The global mesh build is **Python-only**, and the port says so:
  `src/python/eqdyna/MPI4NodalQuant.py:35-41` ("This port keeps the PROVEN
  serial mesh build on every rank and then RESTRICTS it") and
  `driver.py:326-327`.
- **The "~100 s setup" figure has never been measured.** `writeCompTime` is
  declared `= 0` at `globalvar.f90:109`, read once at `eqdyna3d.f90:86`, and
  set to 1 **nowhere** in `src/` or `scripts/` — those two lines are the only
  hits. So `compTimeInSeconds(1)`, the exact quantity the claim is about, has
  never been printed by anything. CLAUDE.md already warns this counter is
  unreachable; this is that warning coming due.

The memory half of the claim survives and is Python's; the time half does not
survive as stated. The cheap discriminating probe (serial `build_solver_state`
at 1 process and at 32 concurrent) costs ~20 minutes and no solver run, and it
should be done before any code is written.

---

## 10. A correction I was given twice and had to refuse twice

Both the handoff and a mid-session instruction told me to "land
`origin/wei/rough-fault-normal-2026-09-22`" for v5.16.0. That branch was
already fully merged before this session began: tip `1d11d89`,
`git log --oneline origin/wei/rough-fault-normal-2026-09-22 ^origin/master` is
EMPTY, and it is an ancestor of origin/master. There was nothing to land;
v5.16.0 is cut from master.

Relatedly, and caught by zofia rather than by me: the handoff credited v5.15.0
with `01a1040`, `d4f625e` and `a16cc98`. **None of the three is in v5.15.0** —
`git merge-base --is-ancestor 01a1040 54e0697` returns non-zero. The published
v5.15.0 notes do not claim them, because I took the release body verbatim from
the release commit's own README block rather than from the handoff's summary.
That is the whole reason the published tag is accurate, and it is the argument
for never writing release notes from a briefing.

All twelve commits credited in v5.16.0's notes were checked the same way:
every one is in HEAD and none is in v5.15.0.

---

## 11. Rule 21: an amendment, not a violation

I was told `87af56e` and `bbf62c8` — two direct commits of mine in the main
checkout — violated rule 21. I read rule 21 in full before accepting that, and
as written it does not say so: *"The main checkout … is the conductor's: no
agent commits in it"*, and *"the conductor merges, agents branch."* Those are
conductor commits in the conductor's own checkout.

The real gap the night exposed is that **rule 21 presumes exactly one
conductor and is silent on two.** Two conductor sessions were live at once and
both had an equally good claim on the main checkout under the rule as written.
Routed to zofia as an amendment, with the note that rule 21's own
`git rev-parse --git-dir --git-common-dir` test cannot tell a conductor from an
agent — so item 68's proposed pre-commit hook would block the legitimate
conductor unless its marker encodes the ownership answer first.

The demonstrated half stands and is the other session's: `88a4012` was a direct
`commit:` into the main checkout, and that same session's 22:12:20 run is the
incident in section 7.

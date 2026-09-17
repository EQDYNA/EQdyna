# Session log — 2026-09-17, v5.10.0 follow-on (wei-lin, autopilot resume)

Resumed from `pathway_forward.md` per this campaign's own rule, not from a
predecessor transcript (predecessor session is no longer addressable; ~780k
tokens of context unrecoverable — expected and accounted for, this is what
the board is for).

## Budget read (stated before spending anything)

24h autopilot window. First hours gated by load, not by clock: a measurement
window (items 33/40) was reported open and explicitly not mine to disturb.
Plan: zero-load work only (reading, scoping, doc landings) until haruto
reports; heavy work (sweeps, oracle re-runs, dispatch that runs `testsys/`)
deferred until then. Tag budget: patch/minor only, direct to master
(project's actual practice — **named deviation from my contract's default of
non-default-branch tagging**, per the brief's own instruction to log rather
than resolve this silently). No major bump, no publish, no force-tag.
Heartbeat while scoping: 600s.

## Liveness check (rule: see before you spawn/assume)

07:54: confirmed haruto alive by PROCESS, not just his relayed message — PID
1591595 (103% CPU) running `eqdyna3d.run_case(..., backend='numpy')` from his
locked worktree `.claude/worktrees/agent-abdc5947b260baa99/`, plus a driver
process fitting per-step scaling (items 33+40 combined: JAX-CPU core scaling
AND the tpv29 numpy-vs-jax gap, same measurement session). 07:58: re-checked,
earlier PIDs had exited but a NEW one was running (n=100 point, "marginal
60->100 s/step" driver) — correctly read as him moving to his script's next
step, not as him being gone; a `ps -p <old pids>` returning empty is not
evidence of absence, a broader `ps aux` grep is. Did not touch
`testsys/perf/`, did not run anything. Load 07:54–07:58: 3.10 -> 2.39 -> 2.92
(his own measurement traffic, not mine).

## Board audit: two discrepancies found between the committed board and reality

1. **Item 29** (stash drop): board row says "resolved 2026-09-16 (drop
   pending owner)". `git stash list` on the current tree is EMPTY — the
   drops already happened, consistent with the launching brief's "Item 29
   CLOSED — I dropped all four dead stashes with owner approval." The board
   row is stale relative to committed reality. Flagging for Zofia to close
   the row on this session's fresh `git stash list` output (empty) rather
   than re-doing anything — nothing to land, the state is already correct.
2. **Item 19(c) / TPV34-35** (P4 row): board text says "specs not even
   fetched ... `scratch/specs/` still has neither spec fetched (step 1, not
   done)". `git log -- testsys/e2e/full_specs.py` shows commit `3677581`
   ("rule 17 step 1: fetch TPV34/TPV35 specs, record full_specs.py entries")
   already on master, part of the v5.10.0 release (HEAD `d47876b`), with
   `NOTES_tpv3435_spec.md` (also tracked, also on master) recording exactly
   what the fetched spec PDFs say. Rule 17 step 1 IS done for both cases;
   the board row is stale. (The raw PDF/txt/html files themselves are
   gitignored per this project's convention like every other `scratch/specs/`
   doc, so a fresh clone's `scratch/specs/` won't show them without
   re-fetching — that is normal, not a regression.) Flagging for Zofia.

Neither discrepancy changed anything I did — re-run-the-evidence caught both
before I acted on the board's stale text, which is the point of the rule.

## TPV34/TPV35 scoping (rule 17 step-3-adjacent, zero-load, no dispatch)

Read-only: grepped `nmat` handling (`meshgen.f90`/`meshgen.py`), grepped for
any CVM/velocity-model ingestion code (none found anywhere in `src/`,
`scripts/`), grepped for external per-node stress-file input (none for a
from-scratch run; only the mode=2 cycle-restart netcdf path exists and it's
a different mechanism), read `swtwNucleation` (`faulting.f90:386-421`).
Findings appended to `NOTES_tpv3435_spec.md` (extended, not a new file):

- TPV35's material need (1D profile per side of fault) does not fit the
  existing depth-only `nmat>1` layering; TPV34's need (real 3D CVM-H) has NO
  existing counterpart at all — the single largest net-new piece either case
  needs. Likely shape: a pre-extracted, decimated, provenance-cited CVM-H
  grid file, matching this project's TPV29 precedent, not a live query lib.
- TPV34's nucleation (shear-stress bump in a circular patch) looks reusable
  from `swtwNucleation`'s existing TPV29/36/37/201/202 family. TPV35's
  nucleation (externally-supplied lowered-yield-stress patch, "built into
  the input data file") has no precedent for a from-scratch run — genuinely
  new.
- No real-recorded-waveform comparison infrastructure exists for TPV35's
  validation angle; per the TPV29 precedent (rule 17 step 6 is a script,
  landed after gating, never itself a gate), this defers and is not a
  blocker for landing the case on cross-code criteria.
- Net read: TPV35 closer on material, further on nucleation/stress-input;
  TPV34 closer on nucleation, further (by a lot) on material. Recommend
  scoping ONE shared material-lookup generalization (fault-relative AND
  external-grid-capable) rather than two bespoke paths, since it touches
  `meshgen.f90`'s material dispatch that every gated case depends on and
  needs its own correctness gate against the 8 existing cases before either
  TPV34 or TPV35 builds on it.

Not dispatched. Still not enough to hand to an implementer blind, but the
next pass has concrete surfaces to read instead of a blank page. Next real
step for TPV35: fetch `tpv35_data_files.zip`. Next real step for TPV34: scope
a CVM-H extraction/decimation approach.

## Deviation note (rule 2 of the loop discipline — code/reality wins, logged not litigated)

Board's own priority order (P1..P4) would put item 39 (TPV30 promotion, P3,
now unblocked) ahead of item 19(c)/42 (TPV34/35, P4). The launching brief
explicitly named TPV34/35 scoping as next. No actual conflict in practice:
TPV30 promotion requires running cases (gating, freezing a reference) —
exactly the heavy-compute class forbidden during the open measurement
window — so it was not actionable regardless of rank. TPV34/35 scoping is
read-only and fits the zero-load window precisely. Recording this so the
rank mismatch doesn't look silently overridden.

## Landings this session

- Docs-only: `NOTES_tpv3435_spec.md` appended (scoping section above).
  This file (session log). No version bump — matches this project's own
  precedent of docs-only commits landing without a tag (e.g. `e845ad0`).
- No code changed. No test tier run (nothing to gate). No tag cut.

## Open at end of this pass

- Waiting on haruto's item 33/40 report — his measurement window is still
  open (confirmed alive by process, not message, as of 07:58). Do not run
  anything under `testsys/` until he reports.
- TPV30 (item 39, P3) is unblocked and is real work once compute is
  available — should be picked up before TPV34/35 implementation per the
  board's own priority, once the load window allows a run.
- Items 35/38 (CI flakes) remain single-occurrence, not investigated
  further per their own rows — no new occurrence found this pass.

## Update: haruto reported, coordinator relayed at ~08:05

Haruto's actual finding on item 33: not just a load problem — TWO real
defects in `testsys/perf/run_numa_scaling.py` would have produced a
misleading answer even on an idle box (F1: generated an all-64-core config,
contradicting the owner's "at most 32, not 64"; F2: no 16-core point,
4x gap between spread-8 and socket-32). He left a corrected driver at
`/tmp/claude-16759/.../scratchpad/item33_driver.py`; promoted its fix
(config-generation only, ceiling/override logic untouched) into
`testsys/perf/run_numa_scaling.py`, verified by `ast.parse` + `--help` (no
idle box yet to run the actual measurement), committed `03ba055`, pushed.
Item 33 itself is STILL not measured — tooling is now trustworthy, the box
is not yet idle enough to trust the number.

Item 40 answered by haruto (relayed, not independently re-run by me):
numpy-pinned 2072 ms/step (linear, 0.2% repeatable), jax-pinned 524 ms/step
(4.0x), jax-unpinned 365 ms/step (5.7x, the board's old "~6x"). Key finding:
unpinned numpy is UNSTABLE (2.08-3.75 s/step, no step-count dependence) while
pinned is linear — CI's 950-1711s tpv29 spread IS this instability
(1711/950=1.80 matches worst-unpinned/pinned ratio almost exactly), so CI's
critical-path number is a placement lottery, not a fixed cost. Optimization
targets identified but NOT applied: `assembleGlobalKU.py:466`
`calcHourglassResist` 38%/step, `:256` `assembleGlobalKU` 29%/step; separately
`meshgen.py:323 build_node_coordinates` is an 8.76s ONE-TIME (not per-step)
cost, backend-independent (identical numpy/jax profiles), flagged as the
cheapest real win precisely because a per-step metric will never show it
moving.

Dispatched `zofia-kaminska` (agent a83d0da3d39a106bb) to write pathway_forward
rows for items 29 (stash drop, done — `git stash list` empty), 19(c) (rule 17
step 1 already landed at `3677581`/v5.10.0, board text stale), 33 (tooling
fixed, measurement still outstanding), 40 (haruto's numbers, marked as
relayed not independently re-verified by me). Not yet returned as of this
edit.

Dispatched `mira-volkov` (agent a0f87ff963afe7929, isolated worktree) to gate
TPV30 (item 39/19(b), P3, unblocked since item 24(b)/(c) landed). Briefed
with exact paths (`scratch/tpv30/case_input_draft/test.tpv30/`,
`case_input/test.tpv29/` as template, `scratch/tpv30/compareTpv29Tpv30.py`
for rule 17 step 6), rule 17's full recipe in order, and a defect found
during my own scoping pass: the old draft plans `par.tpv=36` to reach the
shared nucleation formula (`swtwNucleation`, `faulting.f90:386-421`) instead
of getting its own `TPV==30` branch — the exact impersonation anti-pattern
rule 17 step 3 exists to catch, named explicitly in her brief as a required
fix, not a workaround to carry forward. Not yet returned.

Chose TPV30 over further TPV34/35 work because it is next in the board's own
P3-before-P4 order and does not need any additional data fetch before
dispatch (unlike TPV34/35, which still need `tpv35_data_files.zip` / a
CVM-H extraction approach before an implementation brief could be written
without guessing).

Both agents run in their own worktrees / own file (pathway_forward.md is
Zofia's alone) — no collision: verified by listing what each was briefed to
touch before dispatching a second one.

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

## Update: F3, kai's pruning pass landed, item 19(c) closed (not demoted)

`0b4684c`: F3 landed in `testsys/perf/run_numa_scaling.py` per owner decision
("if the box has ample space, then that's idle") -- replaced the absolute
`load1 <= 1.0` ceiling (refused at load 1.28 on this 64-core box, ~98% idle,
unmeasurable by construction with `train.py` pinned to one core) with
`cpu_busy_fractions()`/`require_idle(cpus, ...)`, gating PER CONFIGURATION on
`/proc/stat` sampled for exactly the cpus about to be used. Refusal without
override now SKIPS just that configuration rather than aborting the whole
run. Verified live on the real contended box, not just `ast.parse`: cpus
[0,1] measured 0-3% busy despite whole-box load 2.4+, correctly passed;
cpu 63 measured 87-100% busy (the actual hot core) and was correctly
REFUSED at the same default ceiling; override still forces a refused config
through. Same commit pins numpy's BLAS threads (`OMP/OPENBLAS/MKL/
NUMEXPR_NUM_THREADS=1`) on all three CI steps that run a python-numpy cell
(item 40 follow-up, one-line env change, no port change).

Ran item 33 for real afterward (my own background process, not an agent):
within-node 1/2/4 cores measured cleanly (616/405/360 ms/step, 1.52x/1.71x
speedup), but within-node-8 and every larger config were correctly SKIPPED
mid-run when specific target cpus went busy from Mira's and Kai's own
concurrent legitimate work -- the fix behaved exactly as designed under real
contention. Verdict on the NUMA-locality question is inconclusive pending a
re-run once the box is quieter; not yet handed to zofia for the board since
it's a partial result and a clean re-run is cheap once available.

`656e66f`: corrected zofia's own immediately-prior uncommitted edit on item
19(c) from "demoted below P4" to "CLOSED, NOT DEFERRED" -- the owner
decision (tpv34/35 sources predate this repo, `scec_archive/` publication
history exists but the implementation does not, `git log --all -S"TPV==34"`
finds nothing) went further than demotion. `SendMessage` was unavailable to
correct her mid-flight, so this landed as a same-day follow-up commit on the
same row rather than a single edit -- named explicitly here so it doesn't
read as two independent events.

`9fe4a54`: kai-fischer's pruning pass (owner-requested, 46:1 added:deleted
since v5.8.2). Comment/blank-line-only pass on `Dockerfile` (71->52 lines),
`.github/workflows/test.yml` (457->428), `.github/workflows/publish.yml`
(126->112). Verified independently before landing (not just his report):
non-comment/non-blank lines byte-identical to current master in all three
files, both before and after cherry-pick; YAML still parses;
`test_ci_workflow_coverage.py`'s own parser still finds the same 21-cell
coverage post-edit. His fresh full sweep (30/30 cells, 1796.5s) was not
re-run here -- two attempts timed out under real concurrent Mira/item-33
load on the shared box -- but the mechanical comment-only proof is the
load-bearing evidence for this class of change and was independently
produced, not taken from his report. Guard audit (task 2 of his brief): all
17 `testsys/regression/*.py` guards reviewed, deletion set EMPTY -- every
guard either has direct catch-evidence or pins a named past incident (rule
8), full verdict table in his completion report. Docs-duplication task
(task 3) deliberately not done by him -- flagged low-confidence rather than
guessed, per his own brief's "if in doubt, leave it."

Mira (tpv30 gate) still in flight as of this update -- confirmed alive by
fresh file activity (`test.reference.results/test.tpv30/frt.canonical.txt`,
a promoted `testsys/parity/evidence_tpv30_vs_tpv29_contrast.py`, and
`src/fortran/faulting.f90`/`src/python/eqdyna/faulting.py` both modified,
consistent with the `TPV==30` branch fix asked for) even though her
`NOTES_tpv30_gate.md` checkpoint is ~35 min stale -- checked broadly with
`find -newer` before treating that staleness as a stall, per the babysitting
rule (file/proc progress, not just the one checkpoint file).

## Open as of this point in the session (superseded further below — see the
## later "Update" sections for what actually happened after haruto reported)

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

## Update: Mira's TPV30 mission — a real STOP-and-report, landed partially

Mira did NOT gate TPV30 (rule 17 steps 4/7 not satisfied) and said so
plainly rather than force a pass. Finding: a full 3-backend sweep showed
python-numpy and python-jax agreeing with each other to ~1e-3 over the whole
3321-node canonical grid at t=20s while BOTH disagree with the Fortran
reference by up to 4.0e8 Pa (~30% relative) at 3310/3321 nodes -- a real,
deterministic port divergence, not per-backend chaos (chaos would not leave
numpy and jax bit-close to each other while both are far from Fortran).
Binary-search localized: bit-exact across all three backends at t=1s (24
steps, confirms item 24(b)/(c)'s new Tv/taper/Drucker-Prager plumbing is
correctly ported), already diverged with a rupture-arrival flip by t=6s
(144 steps). Not narrowed further inside her session budget.

Verified independently before landing (not on her report alone): diffed
`src/fortran/faulting.f90` / `src/python/eqdyna/faulting.py` against current
master -- the only change is the `TPV==30` addition to `swtwNucleation`'s
shared branch list, matching the spec citation (p.16) exactly; grepped every
`case_input/*/user_defined_params.py` and confirmed NO case currently sets
`par.tpv=30`, so this is a true no-op for all 10 currently gated cases.
`testNameList.py`/`testsys/matrix.py` diffs are comment-only (line-by-line
confirmed, `matrix.CASES` unchanged at 10 entries). Ran
`python3 testsys/run.py unit regression` fresh on the merged tree myself:
SUCCESS both tiers. Landed as `a25ed6c`: the branch fix, the new
`case_input/test.tpv30/` compset (unregistered, documents its own gate
status in its README), the Fortran-only reference under
`test.reference.results/test.tpv30/` (not wired to any gate), the promoted
`testsys/parity/evidence_tpv30_vs_tpv29_contrast.py` (report-only), and her
full trace in `NOTES_tpv30_gate.md`. Both worktrees (kai's, mira's) reaped
after confirming clean/merged status; mira's carried a stale lock (PID
2872539, this session's own shared harness process, not an active writer)
and was force-released, recorded here by name and reason.

Item 39/19(b) needs a board update reflecting this (real progress, not
gated, a concrete open divergence with a file:line'd cause) -- routed to
zofia-kaminska next, not written by me.

## Update: item 39/19(b) board row landed; item 33 reframed twice more; owner
## directive on TPV30 validation queued

Zofia's item 39/19(b) row update verified (cited command
`python3 -c "import testNameList; print('test.tpv30' in testNameList.nameList)"`
-> `False`, matches) and landed `8c5bac8`.

**Owner directive, TPV30 (not yet actioned, queued behind the tpv36/37 RSS
run):** Mira's stop-and-report was the right call on the information she had,
but `scec_archive/tpv30/eqdyna-v3.1-{100m,50m,25m}-2015/` (published
submissions, cplot + on-fault stations) already exists and makes rule 17
step 6 satisfiable today. Directive: build the SCEC comparison for tpv30 on
the `evidence_tpv29_scec_comparison.py` pattern (extend, don't duplicate,
mira's own `evidence_tpv30_vs_tpv29_contrast.py`), compute every number from
data (never seed from prose -- tpv29's own Mw 7.45->7.034 correction is the
standing lesson), label which metric-definition rule produced each number. A
500 m agreement is a REGRESSION check, not an accuracy claim (owner's own
tpv30 standing: 100m rank 14/14, 50m 12/14, 25m 8/14 -- coarse ranks last by
construction). If it roughly matches: gate all three backends at 500m with a
measured bound (not THRESHOLD). If not: report it, don't gate -- "a
divergence here would be more valuable than a gate." My own read going in:
her already-found ~30% Fortran-vs-Python divergence by t=6-20s is itself
legitimate grounds to still not gate even if the published-data comparison
looks fine, and the brief needs to leave room for that combined honest
outcome rather than forcing a pass. Then, ONLY if it gates, cost out a 100m
(not 50m) validation run before queuing it -- cheapest resolution with a
real published cross-code comparison point, matching tpv29's own 100m
precedent over its unneeded 50m run.

**Item 33 reframed a second time:** the actual goal is "make jax and numpy
scale well," not produce a locality measurement -- measurement is
instrumental, and I had been routing it at the wrong tool.
`run_numa_scaling.py` (F1/F2/F3, already landed, still worth having) answers
LOCALITY; the scaling curve itself lives in `run_scaling.py`, whose
`PY_THREADS = [1,2,4,8]` structurally stops before 32 while
`FORTRAN_RANKS` goes to 32 -- so the python side was never going to answer
whether either backend reaches 32 cores. Real confound found in that script:
`cores = ','.join(str(c) for c in range(threads))` then bare `taskset` --
on this 8-node x 8-core topology, 16 threads silently spans 2 NUMA nodes,
32 spans 4, and `taskset` binds CPU but not memory (first-touch, may run
fully remote). Directive: extend `PY_THREADS` to match `FORTRAN_RANKS`
(1/2/4/8/16/32, capped at one socket=32 per the owner's own constraint),
replace `taskset` with `numactl --cpunodebind=N --membind=N` (record
placement per point), and run TWO placement policies (compact vs
one-per-node spread) so a knee can be attributed to locality vs
out-of-parallel-work. THEN fix what it finds -- numpy anti-scales (2072ms/step
pinned vs 2075-3747 unpinned, never faster, up to 1.8x slower: thread
oversubscription + page migration is the leading hypothesis, candidates are
explicit BLAS thread control / `--membind` / first-touch init); jax gains
only ~1.44x despite XLA_FLAGS/OMP_NUM_THREADS already being set, so per-step
work size and memory placement are the candidates; `numpy.ufunc.at` (10%,
scatter-add, notoriously serial) is a plausible hard ceiling worth checking
early. Owner routes this to mira-volkov (parity first, then optimize in the
target idiom, re-checking parity every step; all three backends must stay
green at CURRENT bounds -- this is optimization, not a physics change, so no
reference moves). Sequenced explicitly AFTER the tpv36/37 RSS run and the
tpv30 validation, and not to start while the box is loaded -- not yet
dispatched.

**tpv36/tpv37 RSS measurement status:** launched as my own background
process (not an Agent, no async notification), `test.tpv36 x python-numpy`
ran past 14 minutes of CPU time -- not hung (steady ~100% CPU). `sleep`-based
polling is blocked in this environment, so checks were spaced out. Holding
the tpv30 dispatch and the item-33 rescoped dispatch until this finishes,
per explicit instruction that the box cannot take more than one of these at
a time.

**CORRECTED (coordinator caught this before it reached the board): I floated
"the wedge kernel is much more expensive per-step" as the explanation for
tpv36's slowness. WRONG, and retracted -- this is the fifth mechanism this
campaign that did not survive arithmetic, and it was a plausible-sounding
guess offered without doing the arithmetic first, exactly the failure mode
rule 4/6 exist to catch.** The real explanation is pure geometry, no kernel
term needed: tpv36 sets `par.dt = 0.5*par.dz/par.vp` (not `dx`, unlike tpv8),
and `dz = dx*sin(15 deg) = 129.4 m` for its 15-degree dip -- so relative to
tpv8 (dx=500, term=5): timestep is 0.259x -> 3.86x more steps; term 6/5 ->
4.64x steps; z-resolution 3.86x finer -> ~3.86x elements. Predicted work
ratio ~17.9x, observed ~13-14x -- the geometry OVER-explains the slowdown,
so the wedge (`C_degen`) code path is, if anything, slightly CHEAPER per
element-step than the hex path, not more expensive. Recording the wrong
version would have sent item 33's optimization work profiling
`compute_element_shape`'s wedge branch for a cost that isn't there, when the
measured hot spots (`calcHourglassResist` 38%, `assembleGlobalKU` 29%) are
backend-wide, not wedge-specific. Not polling further on the RSS run per
instruction -- the coordinator has a waiter on the process and will resume
this session when it exits.

## Update: RSS measurement re-planned and finished in minutes, not hours

Coordinator caught a second problem before it cost 3 hours: peak RSS is a
setup-time property (mesh/solver arrays preallocate before stepping, do not
grow during it), so a handful of steps reaches the same peak as a full run.
Validated directly rather than assumed: let the in-flight `test.tpv36 x
python-numpy` full run (556 steps, 1415.9s) finish as ground truth (peak
3050252 KiB, also PASSED its bound 7.894064e-09 vs 1e-06 -- confirms it
wasn't stuck, just doing real, geometry-explained work), then ran the SAME
cell truncated to 28 steps (par.term=0.3): peak 3047536 KiB, 0.09%
different. Method validated on the exact cell that mattered, then used for
the remaining three (killed the in-flight full `test.tpv36 x python-jax`
attempt by PID once validated, rather than let ~3h of full-length runs
continue): tpv36-jax 4.34 GB (38.6s), tpv37-numpy 2.91 GB (95.4s), tpv37-jax
4.20 GB (30.1s). All four land in `af892df`: `testsys/matrix.py`
`MEASURED_PEAK_RSS_GB` + `CI_CELLS` (21->25 of 30), `.github/workflows/
test.yml`'s cheap-group job, `CLAUDE.md`'s cell count. Verified before
landing: `test_ci_workflow_coverage.py` PASS (25 cells, no gap/overlap),
fresh `unit regression` SUCCESS.

Proceeding down the queue as planned: TPV30 SCEC-archive validation next
(mira-volkov, dispatching now), item 33's reframed scaling-optimization
mandate after that.

## Update: TPV30 SCEC validation landed (f110c75, c3dace7); item 33 dispatched

Mira's TPV30 validation mission returned cleanly: two separate findings, not
conflated. (A) physics validity ROUGHLY MATCHES the owner's own published
100m 2015 submission (99.2% rupture-extent overlap, area ratio 0.978/1.013,
median slip diff 6.3% over 24 stations, Mw 7.026 computed fresh) -- framed
explicitly as a regression check (the owner's own 100m submission ranks
14/14 among 14 cross-code submissions). (B) port correctness (numpy/jax vs
Fortran, up to 30% divergence by t=20s) UNCHANGED, still unresolved. Gate
action: none, correctly -- both must hold, only (A) does. Independently
re-run by me in the main checkout (scec_archive/ is gitignored, absent from
her worktree; she used a temporary symlink, removed before finishing) before
landing `f110c75` -- numbers reproduced exactly. Worktree reaped after
confirming its uncommitted content matched what I'd already landed
byte-for-byte (not a stale-lock case, just normal cleanup). Board updated by
zofia (`c3dace7`), verified against her own fresh script run before landing.

Dispatched `mira-volkov` (fresh isolated worktree) for item 33's full
reframed mandate: fix `run_scaling.py` (extend `PY_THREADS` to match
`FORTRAN_RANKS` up to 32, replace bare `taskset` with `numactl
--cpunodebind/--membind`, two placement policies -- compact vs spread -- to
attribute any knee, and the SAME per-cpu busy-check discipline as
`run_numa_scaling.py`'s F3 so individual points stay trustworthy even though
this box has not been idle all session), then root-cause and FIX numpy's
anti-scaling and jax's ~1.44x ceiling, checking the two known project
failure modes (HLO literals, jit-inside-a-loop recompilation) first, keeping
all three backends green at CURRENT bounds throughout (optimization, not a
physics change). Not yet returned.

## Update: Mira's item 33 mission returned; a real regression caught by
## gate-axis-3 review before it could land

She delivered honestly and disclosed her own gaps rather than overclaiming:
`run_scaling.py` rewritten (PY_THREADS to 32, numactl replacing taskset --
including a genuinely new numactl-rankfile mechanism for the Fortran/MPI
path since a bare `mpirun --bind-to` can be reissued past an outer numactl
restriction by OpenMPI's own hwloc binder -- smoke-tested, not just written),
compact/spread placement, per-cpu busy-check reused (not reimplemented) from
`run_numa_scaling.py`. Measured (disclosed as taken with `--busy-ceiling 0.9`
override, not the tool's strict 0.2 default, since the box was never
idle): jax gains a real, if modest, 2.35x at 32 cores (cross-checked against
this session's own earlier clean 1-4 core numbers); numpy anti-scales,
mechanism verified by direct microbenchmark (OpenBLAS reads thread count at
`import numpy` time only, confirmed by testing env-var-after-import vs
numactl-pinned-before-start) rather than asserted. jax HLO dumped and read
directly: one `fori_loop`, 1717 fusion + 556 scatter ops on arrays up to
`f64[3818584]` -- real structural parallelism exists, ruling out
"nothing to parallelize" as the explanation for the plateau past 8 cores.
Both known project failure modes (HLO literals, jit-in-a-loop) checked
against the actual code and ruled out, not assumed. Fix landed:
`_narrow_numpy_affinity()` in `eqdyna3d.py`, pins the numpy backend to a
single cpu (bounds the NUMA-migration worst case; does not claim to make
numpy scale, since nothing can with its current single-threaded kernels --
correctly scoped, a real kernel rewrite is out of scope here). New unit
test added. She explicitly flagged that only 1 of 10 numpy cells (tpv8) was
re-verified against its frozen reference before she had to kill her own
`testsys/run.py all` run under this session's box contention (10 concurrent
numpy cells at ~10% cpu each), and recommended completing verification
before landing rather than claiming it done.

**Gate-axis-3 review (wei-lin) caught a real defect her own testing didn't
reach: concurrent numpy cells collide on ONE cpu.** `min(current)` always
picks the SAME lowest-numbered cpu for every process. Verified directly,
not inferred: launched two `python3 -m eqdyna` processes concurrently,
read `/proc/<pid>/status`'s `Cpus_allowed_list` for both -- both `0`. This
exactly explains the ~10%-cpu-each pileup observed during her own 10-cell
parallel sweep, and would have silently made every future concurrent local
sweep (`run_e2e.py --jobs N>1`, this project's own wider gate) serialize
onto one physical core while dozens sit idle -- a real throughput
regression to this project's own verification workflow, not caught by her
brief (which asked her to test the scaling-tool's sequential use case, not
concurrent sweep invocation) or by her one single-cell correctness check
(which doesn't exercise concurrency at all). Fixed directly (small, bounded,
well-diagnosed -- judged not worth a second full dispatch round-trip):
`ordered[os.getpid() % len(ordered)]` instead of `min(current)`, still
narrows to exactly one cpu (NUMA-migration fix unchanged) but different
processes now land on different cpus. Re-verified both ways: her own unit
test rewritten (the old one asserted "always narrows to lowest", which
would have been PID-order-dependent and occasionally wrong under the
corrected behavior -- replaced with a PID-aware assertion plus a new test
simulating two different-parity pids), and the real concurrent-process
check repeated post-fix: two processes now land on cpus 25 and 27, not both
0. All 5 unit tests pass.

Now running (background, `bdqqn40em`) the actual gate: the 9 numpy cells
her mission didn't re-verify (`drv.a6, tpv10, tpv104, tpv1053d, meng2023a,
meng2023cb, tpv29, tpv36, tpv37`) against their frozen references, `--jobs
8` (safe now that concurrent cells no longer collide). Not landing
anything from this mission until that completes green. In progress: 3 of 9
done so far (`meng2023a`, `drv.a6`, `tpv29`), all SUCCESS against their
existing bounds with the fixed affinity code in place.

## Update: item 33 landed (`d70811d`), all 9 remaining numpy cells green

`bdqqn40em` finished: all 9 SUCCESS against their existing frozen references
and bounds, unmoved (`drv.a6` flip-budget 423/450, `tpv29` 9.876230e-15 vs
1e-10, `tpv36`/`tpv37` 7.9e-09/6.9e-09 vs 1e-6, others similarly clean) --
combined with her own `tpv8` check and jax/fortran being provably
unreachable by a change gated on `backend=='numpy'`, this is full coverage.
One more diff/syntax/staleness pass against current master (`d51a2a4`,
zero commits had touched the changed files since her branch point) plus a
fresh `unit regression` before landing.

**Commit-message corruption caught and fixed before it shipped:** the first
commit attempt (`9fd2cd8`) had several backtick-quoted code fragments (`` `mpirun --bind-to` ``,
`` `import numpy` ``, `` `ordered[os.getpid() % len(ordered)]` ``) silently
eaten by tcsh's backtick command-substitution when embedded in a shell
heredoc -- this repo's shell is tcsh, not bash, and backticks are NOT inert
there the way they are in a bash single-quoted heredoc. The commit itself
landed fine (code diff unaffected, message only), but left visible blank
gaps and stray `mpirun`/`import-im6` error noise in the terminal. Rewrote
the message to a file and used `git commit --amend -F <file>` instead of
embedding it inline -- avoids the shell entirely. Re-verified the amended
message contains all three previously-eaten fragments before
cherry-picking. Landed as `d70811d`, worktree reaped clean (nothing
uncommitted, content matched master byte-for-byte before removal).

**Lesson for the rest of this campaign:** never embed a commit message with
backticked code spans directly in a `git commit -m "..."` shell argument on
this box -- write it to a file and use `-F`, every time, regardless of how
short the message looks.

## Update: coordinator found a real spec-level bug candidate in TPV30 --
## separate from the port divergence, dispatched to dunyu-liu

Coordinator fetched `scratch/specs/TPV29_30_Description_v06.pdf` directly
(rule 17 step 1 had never actually been done for tpv30 -- an omission this
whole campaign's tpv30 work carried without anyone catching it until now)
and audited the initial-stress setup line by line. Spec states twice
(lines 28, 424) "the material properties are the only difference between
the two benchmarks" -- same stress tensor, same b-coefficients, same
friction/nucleation. Parameters (Tv, taper, cohesion, bulk friction,
deviatoric ratio, the spec's own 93%-of-yield check) all independently
verified correct. But on-fault initial shear at station
`faultst000dp120` disagrees: our tpv30 measures ~18% HIGH at t=0 (33.05 MPa)
vs both our own tpv29 (27.79 MPa, dx=200) AND the published SCEC v3.1 tpv30
(28.18 MPa, 100m) -- which agree with EACH OTHER. Our tpv30 then visibly
relaxes down to ~28.75 MPa over 2s, consistent with starting ABOVE yield
(spec sets initial state at 93% of yield; an 18% overshoot exceeds it) --
reframing the earlier-documented "relaxation" as a SYMPTOM of a wrong
initial condition, not a resolution artifact. Mechanism candidate (not
established): `par.C_elastic` selects two different on-fault initial-stress
code paths (tpv29 never calls `setPlasticStress`; tpv30 does,
`meshgen.f90:104` region) and the spec requires them to produce an
IDENTICAL tensor.

This directly means, if confirmed: my `d51a2a4` revert was right for a
stronger reason than known at the time (not just "an unexplained
divergence" but a located, spec-verifiable bug candidate), and
`test.drv.a6` (also `C_elastic=0`, currently GATED and SHIPPED) may share
the same code path -- a much bigger finding than tpv30 alone if so.

Dispatched `dunyu-liu` for the controlled experiment (both cases at their
own dx=500 fast-gate resolution, diff the SAME station's t=0 column,
removing the dx=200-vs-100m confound from the coordinator's own
comparison), then code-path localization, then the drv.a6 blast-radius
check -- explicitly barred from fixing this by loosening a bound or
regenerating any reference.

**Rate-limit note, recording per this campaign's own standing lesson:**
first dispatch attempt hit a Fable-5 session-level 429 before even creating
a worktree (nothing to reap). Re-dispatched the identical mission with
`model: "sonnet"` override, matching this repo's own prior precedent
(commit `3677581`'s author note) for the same failure mode.

## Update: dunyu-liu's investigation returned -- CONFIRMED, root-caused, a
## shipped/gated case is affected. ESCALATING, not proceeding further.

Verdict: confirmed, root-caused, NOT fixed. Landed the notes file
(`fda2d4a`) only -- no source change, per the mission's own correct
decision not to fix without an owner call.

**Mechanism, independently verified by wei-lin (not taken from the report):**
`src/fortran/faulting.f90` lines 119/123/127
(`nsdInitTractionVector(1..3)*C_elastic`) and 179/180
(`xyzInitTractionVector(j)*C_elastic`) -- read directly, confirmed present
exactly as reported. For `C_elastic=0` (tpv30, drv.a6), this multiplies the
CORRECT, spec-matching analytic on-fault initial traction by zero, so the
fault's actual initial state instead comes from an off-fault-element-stress
interpolation (`setPlasticStress`'s volumetric tensor, not locally rotated
to the rough fault normal) -- which does not match.

**A real detour that resolved IN FAVOR of the finding, recorded because it
could easily have gone the other way and I want the reasoning on record:**
while independently checking this, `netcdf_io.f90`'s
`netcdf_read_on_fault_eqdyna` appeared to read on-fault init stress from
slot 19/20 while drv.a6's own `case.setup` writes its stated 40e6/-120e6 to
Python-side slots 7/8 -- looked like a SEPARATE writer/reader slot mismatch,
which would have meant dunyu-liu's "same mechanism as tpv30" claim for
drv.a6 was wrong. Traced `scripts/case.setup`'s actual on_fault_vars
netCDF-write remap table (lines 199-227) before concluding anything: it
explicitly writes Python slot 8 -> netCDF `var18` (`init_strike_shear`) and
slot 7 -> `var19` (`init_normal_stress`), which IS what the Fortran reader's
positions 19/20 read. No mismatch -- the analytic value IS set correctly
before the `C_elastic` multiplier discards it. This resolves the detour in
favor of, not against, dunyu-liu's finding.

**Blast radius CONFIRMED directly against the real file, not inferred:**
`test.drv.a6` (`C_elastic=0`, currently GATED and SHIPPED) states 40e6 Pa
shear / -120e6 Pa normal in its own compset
(`case_input/test.drv.a6/user_defined_params.py:134-135`). Its own
committed, currently-PASSING reference
(`test.reference.results/test.drv.a6/faultst000dp075.txt`, t=0.0417s, read
directly) shows **24.6801 MPa shear / 76.8865 MPa normal** -- confirmed
exact match to dunyu-liu's cited figures, from the actual file.

**ESCALATING rather than proceeding.** This is a real, independently
confirmed correctness bug in a case this project has already shipped and
gates on every sweep -- its actual physics does not match its own stated
intent, and the two candidate fixes (drop the `C_elastic` multiplier
entirely, or make `setPlasticStress` locally rotate to the rough-fault
normal) are both real physics decisions that would move `test.drv.a6`'s
currently-passing flip-budget gate and its committed reference. Per rule 7
(reference regen is a deliberate, reviewed commit, never a side effect of a
bug fix) and my own escalation criteria ("a decision changes product
behaviour... rather than just landing a verified fix"), this is not mine
to decide. Not fixing, not regenerating anything, not touching
`test.drv.a6`'s gate. Surfacing this to the coordinator/owner as the
top-priority open item.

## Update: owner decided option B (resolve setPlasticStress onto the local
## fault normal, not drop C_elastic), dispatched to dunyu-liu

Owner independently re-verified all three of my escalated claims (the
`C_elastic` multiplier sites, drv.a6's declared 40e6/-120e6, its actual
24.6801/76.8865 MPa reference) before ruling. Chose option B because the
element tensor is already correct (reproduces the spec's 93%-of-yield check)
and the defect is specifically its RESOLUTION onto the nominal plane instead
of the local rough-fault normal (spec line 436 endorses this architecture).
Flagged that `C_elastic`'s multiplier has no recorded rationale but a
plausible one (avoiding double-counting if the element stress already
reaches fault nodes via the split-node formulation) -- ordered a source
trace of that premise FIRST, before any fix code, with instructions to stop
and report back if it doesn't hold.

Dispatched to `dunyu-liu` with: the Step-1 trace requirement (cheap, must
not be skipped), the local-normal fix using the already-shipped per-node
`dy/dx`/`dy/dz` gradients (`readInputFiles.f90:316`), both backends,
validation against the SPEC (tpv29-vs-tpv30 matched-dx comparison + the 93%
check) rather than against any reference (references currently encode the
bug), explicit prohibition on touching `test.drv.a6`'s reference/gate as
part of this fix (that regen is its own separate, reviewed commit per rule
7, after this lands), and a flag requirement if drv.a6's own compset
declarations become inert for plastic cases (a release-note-worthy
contract change).

**Rate-limit note (third occurrence this session of the identical
pattern):** first dispatch attempt hit the same Fable-5 429 immediately.
Re-dispatched the identical mission with `model: "sonnet"` override.

## Update: Step 1 FALSIFIES option B's premise -- escalating again rather
## than inventing a third fix

Per the mission's own explicit branch ("if the premise is wrong, stop and
report"), dunyu-liu did not write a fix. Verified independently by wei-lin
before accepting this (it overturns an owner-approved plan, so it needed
the same scrutiny as the original escalation): `un`/`us`/`ud` -- the local
rough-fault-normal projection vectors -- are built from `pfx`/`pfz`
(`src/fortran/meshgen.f90:892-905`) gated ONLY by `insertFaultType>0`, with
NO branch on `C_elastic` anywhere in that construction. Read directly,
confirmed. The FE-reconstruction path that projects nodal quantities onto
those vectors (`faulting.f90:89-97`) is likewise unconditional on
`C_elastic`. So the local rough normal is ALREADY applied identically for
both `C_elastic=1` (tpv29) and `C_elastic=0` (tpv30, drv.a6) -- there is no
"resolves onto the nominal plane instead of the local normal" defect for
option B to target. `setPlasticStress`'s element tensor and the analytic
on-fault path were checked term-for-term (by the mission, algebra verified
against this compset's own constants) to already be identical once
resolved onto that same shared local normal, from the same source file.

Landed the notes file only (`72a353b`), no source change -- correctly
matching the mission's own stop condition. Residual mismatch is now
HYPOTHESIZED (not confirmed) as FE-discretization/mass-split reconstruction
error scaling with fault roughness, weakly corroborated by an
already-committed probe (`c7c4f5f`) finding a smaller version of the same
effect at drv.a6's lower roughness -- not chased further, a genuine
numerical-methods design question flagged for routing to a numerics
specialist rather than invented here.

**Not inventing a third fix myself, not re-litigating the owner's decision
in the loop -- surfacing this back to the coordinator/owner as the correct
next step,** same as the original escalation. `test.drv.a6`/`test.tpv30`
gates/registrations/references remain fully untouched.

## Update: owner independently re-verified the falsification, reframed as
## "analytic vs FE-reconstructed initial traction" (a discretization
## question, not a bug), and asked for a 500/200/100m convergence check

Owner traced `meshgen.f90:892` themselves and confirmed no `C_elastic`
branch in 860-930 -- agreed option B had no target. Reframe: with the local
normal applied identically either way, `C_elastic` only chooses WHICH
initial traction the fault sees -- the analytic value (1) or the FE
reconstruction of the element stress tensor (0). The ~18-19% gap is then a
DISCRETIZATION difference (FE reconstruction vs analytic), which should
SHRINK with resolution if that's really all it is. Test: run tpv30 at
500/200/100m, read `faultst000dp120.txt` t=0 shear. Converging closes this
as a documented `C_elastic=0` property; not converging means a real defect.

**Ran it myself (mechanical, no dispatch needed), stated cost first (mesh-gen
dominated once truncated to `par.term=0.2`, ~1-4s at dx=200/dx=100 case.setup
plus a short mpirun since only ~24 steps are needed for the first output
row, not the full 20s duration):**

- dx=500 (already known, dunyu-liu's Step 1): **33.1126 MPa**
- dx=200: **FAILED** -- `case.setup`'s own fault-geometry writer produced a
  self-inconsistent dx=200 decimation (`validateFaultRoughGeometry` rejected
  its own freshly-written `bFault_Rough_Geometry.txt`: derivative columns
  disagree with the surface column by 2.0x the bound at 2 of 202 boundary
  nodes). Reproducible, not a stale-file artifact (timestamps confirm the
  file was written fresh by this exact run). A real, separate finding --
  the geometry decimation path has a bug at dx=200 specifically. Not
  investigated further this pass (out of scope for the physics question);
  worth its own board row.
- dx=100 (native shipped resolution, no decimation, so immune to the dx=200
  bug): **28.0546 MPa**, t=0.00833s first output row.
- Published (SCEC v3.1, 100m, cited by the owner): 28.18 MPa.
- Analytic (bit-identical between tpv29/tpv30's own on_fault_vars_input.nc,
  per dunyu-liu's Step 1): 27.79 MPa.

**33.11 -> 28.05, converging tightly onto the published 28.18 (0.46% off)
and close to the analytic 27.79.** Strong support for "FE-reconstruction
discretization error that shrinks with resolution," not a defect -- though
the middle point (200m) is missing due to the unrelated geometry bug above,
so this is not yet a complete three-point convergence curve. Scratch dirs:
`scratch/tpv30_convergence/test.tpv30.dx100`,
`scratch/tpv30_convergence/test.tpv30.dx200` (failed case.setup, kept as
evidence of the geometry bug).

**Owner directive received mid-run: cut the release first (24 commits sitting
on master), THEN continue this convergence work + the drv.a6 inert-declaration
fix + item-33 requiet-box re-measure, in that order.** Switching to the
release now; the dx=200 geometry bug and completing the convergence
picture are queued right after.

## Update: coordinator reversed the convergence call to CONFIRMED BUG (ran
## the control themselves: tpv29 reads 27.791 at BOTH dx=500 and dx=200,
## resolution-independent since it's the analytic value) -- asked for one
## more localizing run at dx=200 for tpv30, "cheaper than the 500/200/100
## sweep." Dispatched haruto for the release in parallel (not delayed).

Retried tpv30 at dx=200 in a completely fresh directory (ruling out
contamination from the earlier attempt): **identical failure, same numbers**
(`case.setup`'s own `validateFaultRoughGeometry` rejects the freshly-written
dx=200 surface, 2 of 202 boundary nodes over the derivative-consistency
bound by 2.0x, worst at index (0,42), `|d - onesided(y)| = 1.7525e-03` vs
bound `8.9320e-04`). Reproducible, not a fluke.

**Then tried `test.tpv29` at dx=200 myself, to see whether I could reproduce
the coordinator's own reported clean 27.79 MPa control run.** It FAILS
identically -- same error text, same node count, same worst-index, same
numbers, in this environment. I cannot get either case to build at dx=200
right now, so I cannot produce the localizing data point as asked.

**Not accusing the coordinator's number of being wrong -- reporting the
discrepancy plainly and asking, rather than guessing why our environments
disagree.** Possible explanations, not chosen between: a different override
method that bypasses `ensureFaultRoughGeometryForCase`'s writer path
entirely; a different git state; or a genuine intermittent/host-dependent
element in the geometry-validation tolerance check. Either way, this dx=200
geometry-decimation failure is itself real and reproducible in my hands for
BOTH tpv29 and tpv30 -- independent of `C_elastic`, so on its own it neither
confirms nor refutes the initial-stress finding, but it does mean the
"cheaper single run" cannot be obtained here without first resolving (or
working around) this separate validation bug.

Not chasing the geometry bug further without direction -- it's already a
second real, separate defect this session surfaced as a byproduct
(alongside the drv.a6 inert-declaration trap), and digging into it now would
be scope creep on top of an already-large thread. Surfacing it and holding.

Haruto-nakamura (agent `a5be6c04216c5ef56`) dispatched for the minor release
in parallel -- confirmed alive (worktree locked), his own `testsys/run.py
all` gate likely explains the box load spike to 19-31 observed this pass;
not competing with it further.

## Update: TPV30 gate attempted per coordinator instruction, reproduced the
## known divergence, REVERTED rather than forced

Coordinator read commit `f110c75`'s "(A) roughly matches" as the owner's
"sure, gate it" condition satisfied, and asked for the full rule 17
steps 4/7 gate (all three backends, bound from measured worst diff, not
THRESHOLD). Provisionally registered `test.tpv30` (`testNameList.py`,
`testsys/matrix.py`, bound `1e-6` explicitly marked PROVISIONAL, `GATE=
'abs-max'`) and ran a fresh 3-backend sweep against the existing frozen
reference to find out, rather than assume either outcome.

**Result: the known divergence reproduced exactly.** `test.tpv30 x fortran`
SUCCESS (bit-exact, `max|diff|=0.0`, unaffected as expected).
`test.tpv30 x python-numpy` and `x python-jax` both FAILED identically:
`max|diff|=4.035089e+08` (same row/col, same magnitude as Mira's original
finding -- this is the SAME bug, not a new one). Per rule 2/5 and my own
stated commitment before running it, **did not gate** -- reverted
`testNameList.py`/`testsys/matrix.py` to their committed state (`git
checkout --`), confirmed clean (`test.tpv30 in nameList` -> `False` again).
The coordinator's instruction was reasonable given what commit `f110c75`
said in isolation, but condition (B) (port correctness) was never resolved
by the SCEC validation and this fresh run confirms it directly -- (A)
matching does not make (B) go away, exactly the distinction the row was
written to preserve. Nothing to land from this attempt; the outcome IS
the value (a second, independent confirmation of the divergence, on a
fresh run, is stronger evidence than the first one alone).

**Item 33 mandate escalated again:** not a curve to report -- an
optimization to land. Owner: "if not, go optimize it" / "at least, jax
should scale really well." jax's measured 1.44x (1 core to unrestricted) is
now treated as a DEFECT to root-cause, not a property to report, since XLA-CPU
parallelizing badly when configured and fed properly is itself anomalous.
Candidate causes, cheapest first, none to be believed until reproduced: (1)
memory placement -- `numactl --cpunodebind=N --membind=N` over bare
`taskset` (leading hypothesis, both backends); (2) latency- vs
throughput-bound -- jax already wins ~4x at one core via FUSION, so if the
step is a long dependent-op chain there is little intra-op parallelism to
exploit regardless of core count; dump the HLO (op count, tensor shapes)
BEFORE optimizing anything, since this would reframe the fix as exposing
parallelism (batching/vmap/scan) rather than tuning threads; (3) scatter-add
-- `numpy.ufunc.at` at 10%/31 calls per step, whose jax equivalent
(`.at[].add()` on duplicate indices) is the same op this repo's own notes
say lowers to atomics on GPU and may serialize outright on CPU, a plausible
hard scaling ceiling independent of core count; (4) HLO literals -- this
repo's own documented lesson, a jitted function closing over an array pays
it as an HLO literal (~7x); check jitted entry points pass arrays as
ARGUMENTS; (5) recompilation -- also a documented lesson, confirm `jax.jit`
is built once, never inside a run function. Items 4/5 are known project
failure modes, check first. Numpy target: stop ANTI-scaling (currently never
faster than 1 core, up to 1.8x slower) -- `ufunc.at` + BLAS thread
oversubscription are the leads. Mandate: mira-volkov measures, root-causes,
FIXES, re-measures, keeps all three backends green at CURRENT bounds (no
reference moves, no bound loosens) -- if jax is genuinely at a real ceiling
for this problem shape that's acceptable, but it needs HLO evidence, not an
assertion. Still queued behind tpv36/37 RSS + tpv30 validation, not yet
dispatched.

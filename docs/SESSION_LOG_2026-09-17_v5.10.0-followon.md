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
   `docs/notes/NOTES_tpv3435_spec.md` (also tracked, also on master) recording exactly
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
Findings appended to `docs/notes/NOTES_tpv3435_spec.md` (extended, not a new file):

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

- Docs-only: `docs/notes/NOTES_tpv3435_spec.md` appended (scoping section above).
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
`docs/notes/NOTES_tpv30_gate.md` checkpoint is ~35 min stale -- checked broadly with
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

Haruto's item 33/40 findings landed as commit `03ba055` (two tooling defects
in `run_numa_scaling.py` fixed; item 33 itself still not measured; item 40's
numpy/jax per-step numbers and optimization targets identified) -- full
figures in pathway_forward.md items 33 and 40 (relayed, not independently
re-run by me, flagged as such there).

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
plainly rather than force a pass. Finding (full numbers in
pathway_forward.md item 19(b)/39): a real, deterministic port divergence
between the Fortran reference and both Python backends, binary-search
localized to between t=1s (bit-exact all three) and t=6s (already
diverged). Not narrowed further inside her session budget.

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
full trace in `docs/notes/NOTES_tpv30_gate.md`. Both worktrees (kai's, mira's) reaped
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
conflated -- (A) physics validity roughly matches the owner's published 100m
2015 submission, (B) port correctness (numpy/jax vs Fortran) unchanged,
still unresolved; full numbers in pathway_forward.md item 39's `f110c75`
update. Gate action: none, correctly -- both must hold, only (A) does.
Independently re-run by me in the main checkout (scec_archive/ is gitignored, absent from
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
`run_scaling.py` rewritten (PY_THREADS to 32, numactl replacing taskset,
compact/spread placement, per-cpu busy-check reused from
`run_numa_scaling.py`) and `_narrow_numpy_affinity()` landed in
`eqdyna3d.py`. Measured numbers (jax 2.35x at 32 cores, numpy anti-scaling
mechanism, HLO dump) are in pathway_forward.md item 33's own update for this
commit, not repeated here. She explicitly flagged that only 1 of 10 numpy
cells (tpv8) was re-verified against its frozen reference before she had to
kill her own `testsys/run.py all` run under this session's box contention,
and recommended completing verification before landing rather than claiming
it done.

**Gate-axis-3 review (wei-lin) caught a real defect her own testing didn't
reach: concurrent numpy cells collide on ONE cpu** (`min(current)` picked the
same lowest cpu for every process -- full repro/fix detail in
pathway_forward.md item 33). Fixed directly (small, bounded, well-diagnosed
-- judged not worth a second full dispatch round-trip) rather than sent
back. Re-verified both ways (rewritten unit test, repeated concurrent-process
check). All 5 unit tests pass.

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
and audited the initial-stress setup line by line against it. Finding
(numbers and mechanism candidate in pathway_forward.md item 39's row): our
tpv30's on-fault initial shear at `faultst000dp120` disagrees with both our
own tpv29 and the published SCEC v3.1 tpv30, which agree with each other,
then visibly relaxes toward the target over ~2s -- reframing the
earlier-documented "relaxation" as a symptom of a wrong initial condition,
not a resolution artifact.

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
confirmed directly by reading `src/fortran/faulting.f90` (the
`*C_elastic` multiplier zeroing the analytic on-fault initial traction for
`C_elastic=0` cases) -- full mechanism and the drv.a6 declared-vs-actual
numbers are in pathway_forward.md item 39's "the stakes" paragraph, not
repeated here.

**A real detour that resolved IN FAVOR of the finding, recorded because it
could easily have gone the other way and I want the reasoning on record:**
while independently checking this, a suspected netCDF slot mismatch between
`case.setup`'s writer and `netcdf_io.f90`'s reader for drv.a6's on-fault
init stress looked like it might invalidate dunyu-liu's "same mechanism as
tpv30" claim for drv.a6. Traced `scripts/case.setup`'s actual remap table
(lines 199-227) before concluding anything: no mismatch -- the analytic
value IS set correctly before the `C_elastic` multiplier discards it. This
resolves the detour in favor of, not against, dunyu-liu's finding.

**Blast radius CONFIRMED directly against the real file, not inferred:**
`test.drv.a6` (`C_elastic=0`, currently GATED and SHIPPED) is affected --
its declared compset values and its actual passing reference disagree,
matching dunyu-liu's cited figures exactly (see item 39's row for both
numbers).

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
the same scrutiny as the original escalation): the local rough-fault-normal
projection is already applied identically for `C_elastic=1` and `C_elastic=0`
cases -- no branch on `C_elastic` exists anywhere in that construction (read
directly, confirmed; see pathway_forward.md item 39's exonerated-hypothesis
list for the file:line detail) -- so there is no "resolves onto the nominal
plane instead of the local normal" defect for option B to target.

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
row, not the full 20s duration):** dx=500 (already known from dunyu-liu's
Step 1), dx=200 FAILED (a separate, reproducible `validateFaultRoughGeometry`
rejection of its own freshly-written dx=200 decimation -- same bug pinned in
pathway_forward.md item 39's row, affects both tpv29 and tpv30 identically),
dx=100 succeeded. **The dx=500->dx=100 trend converges tightly onto the
published value** (full numbers in item 39) -- strong support for
"FE-reconstruction discretization error that shrinks with resolution," not a
defect -- though the middle point (200m) is missing due to the geometry bug
above, so this is not yet a complete three-point convergence curve. Scratch
dirs: `scratch/tpv30_convergence/test.tpv30.dx100`,
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

## Update: release attempt #1 hit a session rate limit mid-sweep,
## re-dispatched with model override; coordinator now running the dx=200
## localizing measurement themselves, correctly not duplicated here

`a5be6c04216c5ef56` failed on a Fable-5 session-limit 429 -- per its own
final note, it was "well into the e2e sweep with real successes" when cut
off. Worktree was gone by the time I checked (auto-cleaned on termination),
so nothing to resume from or reap -- re-dispatched fresh, not continued.
**Fourth rate-limit hit this campaign, same pattern each time**: re-dispatched
identical mission with `model: "sonnet"` override (agent `a92c110247686a6e9`),
recorded here per the coordinator's own standing instruction to record which
model ran what.

Coordinator is running the dx=200/dx=500 `test.tpv30` localizing measurement
(t=0 initial shear) themselves and will relay the number -- explicitly told
me not to duplicate it, so the geometry-decimation-bug thread and my own
dx=200 attempts are paused here, not abandoned. `case.setup`'s
inert-on-fault-stress-declaration warning (drv.a6's now-confirmed-live trap)
remains queued behind the release landing.

## Update: release attempt #2 (`a92c110247686a6e9`, sonnet) completed its
## audit/gate but correctly stopped before committing -- assembled the
## rest myself rather than a third full dispatch

Full green gate confirmed (fresh, on pristine pre-release master `bd9d81e`):
`./install-eqdyna.sh -m ubuntu` clean, `python3 testsys/run.py all` --
unit/regression SUCCESS, **e2e 30/30 cells SUCCESS**, wall clock 2263.3s.
CI independently confirmed green on `72a353b` (the last code-touching
commit; the 3 commits after it are docs-only, correctly excluded from CI by
`paths-ignore`).

**Caught and corrected three errors in MY OWN dispatch brief before they
reached the release notes** -- verified independently by wei-lin, not taken
on the report alone: (1) item 39's bound tightening already shipped in
v5.10.0 (`git show v5.10.0:testsys/matrix.py` confirms the "TIGHTENED
2026-09-17" comment predates this release) -- excluded from the v5.11.0
note. (2) The Docker.guide.md de-pinning commit (`ca83335`) is also an
ancestor of v5.10.0 (`git merge-base --is-ancestor` confirmed) -- excluded.
(3) `scratch/specs/` is gitignored repo-wide (`.gitignore:3`), so nothing
there is ever "committed" -- corrected to cite the spec by name/part in
already-committed prose (`docs/notes/NOTES_tpv30_gate.md`, `testsys/e2e/full_specs.py`)
instead.

**Correctly stopped before committing anything:** rule 15 requires the
`pathway_forward.md` Tasks-done row in the SAME commit as the VERSION bump;
since `pathway_forward.md` is zofia's file, committing without her row first
and then tagging would have hit `test_release_complete.py`'s
`check_pathway_tasks_done_row` guard on the TAG-triggered CI run --
a guaranteed, self-inflicted red, correctly refused rather than risked.

His worktree turned out to hold no uncommitted file changes (the VERSION/
README content was drafted in his report text, not written to disk) --
rather than a third full dispatch just to write already-agreed content,
assembled it myself: dispatched `zofia-kaminska` for the Tasks-done row
(landed, diff verified: exactly one row, newest-first order preserved, no
other file touched), then wrote `README.md`'s new v5.11.0 news block
(terse, user-facing, per rule 15's own style guidance -- full detail stays
in `pathway_forward.md`/the GitHub Release body), moved v5.10.0's block
into `pastReleaseNotes.md`, bumped `VERSION` and the Fortran banner string
in `src/fortran/eqdyna3d.f90` to match. Rebuilt clean, full gate
(`python3 testsys/run.py all`) running now in the background before the
release commit.

## Update: coordinator fully explained the tpv30 mechanism -- NOT a bug,
## a documented property of C_elastic=0. Corrected the release wording
## before it shipped.

Full explanation, verified by the owner via source read (mechanism and every
supporting number now in pathway_forward.md item 39's row): elastic cases
see the exact analytic projection, plastic cases see an FE reconstruction of
`setPlasticStress`'s field instead (the analytic value is computed and
correctly discarded to avoid double-counting), and the ~19% gap is two
valid representations of the same physical state on a rough mesh, not a
defect -- this explained every observation six earlier hypotheses hadn't.
Fault morph also verified exact -- geometry ruled out as an error source.

**Corrected before shipping, not after:** edited `README.md`'s v5.11.0 news
block myself (still uncommitted at that point) from "confirmed defect" to
"documented property of C_elastic=0, not a bug -- natural fix is a finer
gate resolution, not a code change." Dispatched `zofia-kaminska` (agent
`a953ba2d6e0652ce5`) to correct the SAME characterization in
`pathway_forward.md` (both the freshly-added v5.11.0 Tasks-done row and
item 39/19(b)'s own row), fold in two new findings (dx=200 blocked for
BOTH tpv29 and tpv30 identically -- not tpv30-specific, a legitimate
resolution refused by a marginal edge-tolerance issue in
`validateFaultRoughGeometry`; and a non-urgent design option for an EXACT
tpv30 initial traction via `T_init = analytic - FE reconstruction`, owner's
call, not landed), and explicitly preserve the drv.a6 inert-declaration
finding unchanged (that one was never in doubt). Not yet returned.

Release gate (`br2ekcqll`) still running in the background, now past
regression into e2e.

## Update: owner ruled "no need for finer grids" -- gate-it-finer wording
## retracted; a real fix direction (discrete-equilibrium correction) given
## and dispatched to dunyu-liu

Owner: "gate it finer" is REJECTED, not just softened -- tpv29 is itself
FE-reconstructed at dx=500 after t=0 and shows no comparable error there, so
mesh coarseness cannot explain the initial state specifically. New, sharper
framing (full reasoning in pathway_forward.md item 39): a uniform far-field
tensor on a faceted rough fault is not in discrete equilibrium, so FE
reconstruction disagrees with the exact analytic projection `case.setup`
already computes and discards. Fix direction (not a specified fix, needs a
numerics read first): carry `(analytic projection - FE reconstruction at
t=0)` as a correction offset instead of zeroing it.

Corrected `README.md`'s v5.11.0 note myself (removed "gate finer"/dx=100
wording, replaced with the mechanism + "fix is a design question in
progress" framing, still uncommitted pending the release). Dispatched
`zofia-kaminska` (agent `a89c1013cdabc87f5`) for the matching
`pathway_forward.md` correction (surgical -- remove finer-grid conclusion,
add the discrete-equilibrium direction and the exact success criterion:
tpv30 matches tpv29's 27.791 MPa AT dx=500, no pre-arrival relaxation).
Dispatched `dunyu-liu` for the actual design+implementation, spec in hand
(`scratch/specs/TPV29_30_Description_v06.pdf`) -- first attempt hit the
same Fable-5 429 immediately (fifth occurrence this campaign, same
pattern), re-dispatched with `model:"sonnet"` override (agent
`a2d9003fa3ce1790f`). Both not yet returned.

## Update: owner retracted the discrete-equilibrium framing too -- it IS a
## confirmed bug, tightly localized to the strike-direction traction
## extraction, ratio 1.18926. Routing corrected from dunyu-liu to
## lars-eriksson (mechanical bug hunt, not new-physics research).

Sharper evidence this round (full numbers in pathway_forward.md item 39's
"established by measurement" list): every input checked identical between
tpv29 and tpv30 except the shear component, off by exactly ratio 1.18926 --
isolating the defect to `us` or the strike-projection of nodal forces
specifically. Seven hypotheses already measured and rejected.

**Coordinator corrected the routing mid-thread**: this is `lars-eriksson`'s
surface (bounded expression, known-correct target, sign/convention error to
locate at file:line, no fix proposed by contract) not `dunyu-liu`'s
(research-heavy new implementations with no reference). Dispatched
`lars-eriksson` (agent `a76f4c9df9083a1f6`) with the full evidence chain and
the explicit target ratio to explain (1.18926) -- read-only, no collision
with the still-in-flight `dunyu-liu` mission (`a2d9003fa3ce1790f`, worktree
confirmed still locked/alive), which I am letting run to completion rather
than guess whether to cancel: if lars locates the real bug first, dunyu-liu's
discrete-equilibrium-correction design becomes a workaround for a problem
that has a direct fix instead, and I'll evaluate both reports against each
other and against the target number when they return, not act on either
alone.

Corrected `README.md` a second time (still uncommitted, pending the
release) and dispatched `zofia-kaminska` (agent `ad7523e125fe386e7`) for a
third `pathway_forward.md` correction on the same thread -- "not a defect"
retracted, "confirmed bug, ratio 1.18926, localization dispatched to
lars-eriksson" is now the standing wording. Release gate (`br2ekcqll`)
down to single-digit remaining cells, still green so far.

## Update: lars-eriksson's localization was WRONG -- caught it myself before
## relaying it, then the coordinator found the REAL root cause independently

`lars-eriksson` reported the bug was in `us`'s normalization
(`meshgen.f90:896-897`, missing a `pfz` term), predicting the measured
1.18926 ratio. Checked his arithmetic myself before relaying it: `us` as
ACTUALLY shipped is self-normalizing, not the vector his prediction
requires -- his number came from a hypothetical alternate formula he never
confirmed was the shipped one. Decisively: the analytic Python path uses the
IDENTICAL construction and produces tpv29's confirmed-correct value -- if
this formula were the bug, tpv29 would be wrong too, and it isn't (matches
the disproof recorded in pathway_forward.md item 39's exonerated list). Did
not relay this as the located bug; reported the disproof back instead.

**Coordinator then found the REAL root cause independently: tpv30 is
running a half-implemented SPEC METHOD, not a numerics bug in a shared
formula.** `TPV29_30_Description_v06.pdf` Part 8 offers two mutually
exclusive initial-condition methods. Method 1 (stress CHANGE from initial,
no gravity, MANUALLY SPECIFIED fault traction, no boundary tractions
needed) is what tpv29 runs. Method 2 (TOTAL stress, explicit gravity,
fault traction DERIVED from the field, REQUIRES boundary tractions if the
mesh boundary can move) is what tpv30 runs -- confirmed by wei-lin
independently: the `(1-C_elastic)` gravity body force at
`assembleGlobalKU.f90:20` exists exactly as described, and a grep for
`boundaryForce`/`boundary_force`/`BoundaryTraction` across both Fortran and
Python returns NOTHING -- Method 2's required boundary tractions are
genuinely unimplemented. Without them the initial state is not in static
equilibrium, the medium deforms from t=0, and the fault traction is read
off a body already in motion -- exactly the measured 33.05->28.75
pre-arrival relaxation.

**Owner decision: move tpv30 to Method 1** (not "complete Method 2") --
drop the gravity body force for this case, keep `setPlasticStress`'s
element field for the yield calculation ONLY, and use the manually-specified
on-fault initial traction `case.setup` already computes and writes
(`on_fault_vars[7]/[8]/[49]`) but currently discards via the `C_elastic`
multiplier. This retires the earlier double-counting objection: under
Method 1 the element stress never generates fault traction at all, so
applying `T_init` directly is the method working as specified.

`dunyu-liu`'s currently-in-flight mission (`a2d9003fa3ce1790f`) was briefed
on the now-superseded "discrete-equilibrium correction" framing
(`T_init = analytic - FE reconstruction`) rather than this cleaner
"switch to Method 1, apply T_init directly" fix. Cannot redirect mid-flight
(no SendMessage) -- will evaluate his report against the sharper Method-1
spec when it returns rather than guess now; a corrected follow-up dispatch
is likely needed regardless of what he produces. Release gate still green,
8 processes remaining.

## Update: release gate went green (30/30), release committed and pushed;
## dunyu-liu's offset mission returned -- verified numerically correct but
## architecturally wrong, sent back for the literal Method 1 switch

Local gate (`br2ekcqll`) finished: **30/30 e2e cells SUCCESS**, unit+
regression SUCCESS, wall clock 2576.2s. Found and fixed one more stale spot
before finalizing: README.md's tpv30 bullet still described the disproven
lars-eriksson "strike-traction" framing -- corrected to the actual
half-implemented-Method-2 explanation before committing. Assembled and
pushed the release commit (`2e42b4b`, VERSION 5.10.0->5.11.0, README news
block, pastReleaseNotes.md move, Fortran banner string) after confirming no
commit had touched those 4 files since I started editing them (zero
staleness). CI running on `2e42b4b`; fortran-a/b, unit-regression, and build
already green as of this checkpoint, the three python e2e jobs (meng,
cheap, tpv29) still in progress -- holding the tag until all are green.

`dunyu-liu`'s discrete-equilibrium-correction mission (`a2d9003fa3ce1790f`)
returned. Independently verified myself, not just his report: his fix is a
persistent correction-offset (computed once at step 1, added back every
step) with gravity left on, not the decided Method-1 architecture -- stable
and matching the target numerically, but confirmed via diff that
`assembleGlobalKU.f90`'s gravity gate is byte-identical to before. Reported
this plainly rather than landing a numerically-correct-but-differently-built
fix or silently redoing it myself -- this is a "what does the physics
represent" decision, not mine to make. Owner ruled: **do the literal Method
1 switch**, backed by new evidence that the FE force assembly itself is
correct (full re-derivation numbers in pathway_forward.md item 39's
"established by measurement" list) -- 33.05 MPa and 27.79 MPa are two
different physically-valid quantities that only converge as dx->0, not a
bug in the assembly. Re-dispatched
`dunyu-liu` (agent `a50faafb5e0d2f44e`) fresh for the literal switch: drop
the gravity body force for `C_elastic=0`, keep `setPlasticStress` for yield
only, apply `T_init` directly with the `*C_elastic` gating REMOVED (not
offset-corrected), both backends, verified against the same full-window
stability check I already validated, plus drv.a6's declared 40/-120 MPa
confirmed to go LIVE under Method 1 (a real, expected behavior change).
Explicitly told to leave the rejected offset worktree untouched -- I reap
it only once the Method-1 build is verified, not before.

## Update: v5.11.0 TAGGED AND RELEASED. dunyu-liu's Method-1 build (also
## rejected per the "leave it open" decision) returned working but is not
## landed; one incidental, real, currently-harmless finding extracted and
## recorded separately. Both tpv30 worktrees reaped.

CI went green on `b63b959` (fresh, own run). Confirmed every commit added
after it (`e731c25`, `d74e180`, `9999928`, `7028fc7`, `e588708`) touches
only `docs/SESSION_LOG_*.md` and `pathway_forward.md` -- both provably
CI-path-ignored and behavior-inert -- so `e588708`'s code state is
byte-identical to the CI-validated `b63b959`; tagged `e588708`, not the
stale `b63b959`. Tag + `gh release create` run as one uninterrupted chained
action per rule 15's amendment (`v5.11.0` -> `e588708`). Two same-second
duplicate push-triggered CI runs fired on the tag push -- matches item 35's
already-documented, not-investigated-further pattern exactly, not a new
concern.

`dunyu-liu`'s Method-1 mission (`a746fa7807cf3485e`) returned AFTER the
"leave it open" decision -- reported a working, fully-verified
implementation (27.79127 MPa matching tpv29 bit-identically, drv.a6 reading
its own declared 40/-120 MPa, 27/30 sweep green with only drv.a6's 3 cells
failing as anticipated). **Not landed, per the owner's final decision** --
the whole approach was retracted before this report even arrived, based on
newer evidence (the resolution-independent published data) that this
report's own author never saw. Worktree reaped after extracting one
genuinely separate, valuable finding: his own report flagged an incidental
numpy aliasing bug at `driver.py:182` (`xp.asarray(inv['stress_i0'])`
missing a defensive `.copy()`) that broke `test.tpv8 x python-numpy` by
1.2e17x in HIS build. Verified myself before accepting or extracting
anything: in CURRENT MASTER, `inv['stress_i0']` is read exactly ONCE
(this same line) and nothing else re-reads it -- so the missing copy is
REAL but CURRENTLY HARMLESS; it only became active because his own
(rejected) code added a second read later in the step function that his
own work needed. Not landing the one-line defensive fix myself either --
it emerged entirely from rejected work and the owner's instruction was "do
not dispatch a fix," not "extract salvageable pieces" -- dispatched
`zofia-kaminska` for a new, separate, low-priority pathway_forward.md row
instead, so it isn't lost for whenever someone next needs a second read of
that value. Both tpv30 investigation worktrees (`a2d9003fa3ce1790f` offset,
`a746fa7807cf3485e` Method-1) now reaped.

Sixth occurrence of the identical Fable-5 429 pattern this campaign, on the
very first attempt at this mission. Re-dispatched immediately with
`model:"sonnet"` override (agent `a746fa7807cf3485e`), no orphaned worktree
left behind by the failed attempt. CI on the release commit (`2e42b4b`, run
`35279131717`) reached all-green (build, both fortran groups,
unit-regression, and all three python jobs) before the next reversal below
made its own release-note wording stale.

## Update: the owner's own published cross-resolution data RETRACTS Method
## 1 too -- Method 2 is proven correct, this is a genuine fault-node
## force-balance defect, not an architecture question

Coordinator compared the owner's three published TPV30 submissions
(100/50/25m): resolution-independent, showing NO relaxation at any
published resolution (full numbers in pathway_forward.md item 39's row).
This build's 500m run relaxes visibly instead. Two kills: a discretization
artifact would grow with element size (it doesn't); a discretization
artifact is a static offset, not something that relaxes toward the right
answer over time (this one does) -- relaxation means the system isn't
starting in equilibrium. Method 2 (what tpv30 already, correctly, runs per
spec) demonstrably WORKS and produces a properly equilibrated state at every
resolution the owner published. Switching to Method 1 would abandon a
proven method to dodge a defect nobody has actually located yet.

**Both prior fix attempts retracted before landing, worktrees KEPT (not
reaped) per instruction:** the discrete-equilibrium correction-offset
(`agent-a2d9003fa3ce1790f`) and the literal Method-1 switch
(`agent-a746fa7807cf3485e`, still in flight when this landed -- will not
act on its report when it returns). Corrected `README.md`'s v5.11.0 note a
THIRD time -- since `2e42b4b` was already pushed and CI was tracking it,
this required a genuine follow-up commit (`b63b959`), not an amend; CI
re-triggered fresh on it (README.md is deliberately in the CI path-trigger
set, unlike `pathway_forward.md`/`PROJECT_RULES.md`) -- holding the tag for
THIS commit's green, discarding `2e42b4b`'s now-stale green as the tag
target. Wrote the new wording deliberately conservatively (states the
established facts -- resolution-independence, no relaxation in the
published data, everything upstream of the fault nodes exonerated -- without
committing to an unconfirmed specific mechanism) given how many times this
exact paragraph has had to be rewritten today.

Dispatched, in parallel: `zofia-kaminska` (agent `abda203611685ab45`, edits
master directly, no isolated worktree) for the pathway_forward.md
correction (sixth pass, retracting Method 1/2, recording the fault-node
force-balance reframing and the class-split residual data); `lars-eriksson`
(agent `ac82ac98e71d89979`, isolated worktree) for the NEW, sharper
localization task -- fault-node force imbalance in the C_elastic=0 path,
using the class-split residual probe (interior max 0.0066, boundary max 11,
near-fault max 594/mean 5-13.5) as the lead, explicitly barred from
re-chasing the now-exonerated list (tensor, gradients, morph, arn/masses,
un/us/ud, the strike-vector formula, the FE assembly algorithm itself).
Both not yet returned.

## Update: lars-eriksson's SECOND localization attempt also wrong -- caught
## it myself again, with a cleaner argument this time

He reported `replaceSlaveWithMasterNode` (`meshgen.f90:760-773`) as the
defect. Rejected, two independent reasons, neither requiring a re-run (both
now recorded in pathway_forward.md item 39's exonerated-hypothesis list):
(1) the check is geometrically correct for this fault's nominal y=0 geometry;
(2) the routine has no `C_elastic` dependency, so any defect here would
already show up in the confirmed-correct tpv29 too, and it doesn't -- also
corroborated by the bit-identical mass measurement already on record. Not
relayed. Root cause remains unlocated -- four plausible mechanisms have now
failed against direct measurement. CI on the release commit
(`b63b959`) still in progress; `dunyu-liu`'s retracted Method-1 mission
(`a746fa7807cf3485e`) still in flight, will not be landed when it returns.

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

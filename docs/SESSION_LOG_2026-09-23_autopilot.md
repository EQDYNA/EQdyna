# Session log — 2026-09-23 autopilot (wei-lin)

Resumed 2026-09-22 23:43 from a clean handoff: origin/master `894cdc1` = tag
`v5.16.0`, VERSION 5.16.0, clean tree, three worktrees, box at 12/64. All three
releases of the previous session were done; nothing was re-cut.

Autonomous grant, stated back as required: minor and patch tags on `master` of
THIS repo, merge and tag authority on `main`/`master` directly. No major bump,
no force-update or rewrite of an existing tag, no publish to a package index,
nothing outward-facing beyond this repo.

## Headline

| | |
|---|---|
| landed | **v5.16.1**, `5b7a278` — pathway items **70** and **68**, both mechanical enforcement, no solver change |
| measured | item **73**'s 16/32-rank remainder — **the non-blocking halo is declined at every rank count that exists** |
| measured | item **64** Stage 0 — `build_solver_state` is **4.42 s** serial, **7.92 s** mean at 32-way concurrency, **1.65 GB per rank** |
| refuted | the inherited "~100 s setup" figure, and the wave argument above 12 ranks |
| NOT started | Stages 1-3 of the setup rewrite — reporting first, as instructed |

## 1. Landing — v5.16.1 (`5b7a278`), 2026-09-23 00:14

Two missions, dispatched concurrently into their own worktrees, merged one at a
time behind my own gate.

**Item 70 — the e2e run-tree lock.** `testsys/runlock.py`, an exclusive
non-blocking `flock` taken at Gate 0 of `run_e2e.py:558` (ahead of the build,
not merely ahead of the rotation) and `run_e2e_full.py:113`, plus
`testsys/regression/test_e2e_run_tree_lock.py` (9 checks, mutation-verified by
the mission). This is the mechanism behind rule 21a and it closes the failure
that cost the v5.16.0 gate 40 minutes: a second sweep renaming the first's live
`test/` tree out from under it, producing a `FileNotFoundError` at 1500.1 s that
printed `FAIL` and read exactly like a parity regression.

My own verification, run fresh rather than taken from the report — hold the lock,
then invoke the REAL entry point:

    holder pid 2715746
    second invocation exit 1 in 0.2s
    ...
    NOT waiting. NOT rotating anyway. NOT falling back to a different directory --
    each of those is the silent fallback rule 2 forbids.
    BUILD ATTEMPTED? False

0.2 s and no build: the refusal happens before anything is spent.

**Item 68 — the main-checkout pre-commit hook.** `testsys/hooks/pre-commit`
refuses when `--git-dir` equals `--git-common-dir` (both resolved with `pwd -P`
first, because git returns the bare string `.git` for both at the top of the main
checkout); `install-eqdyna.sh` sets `core.hooksPath testsys/hooks`; guard at
`testsys/regression/test_precommit_main_checkout_guard.py`. Equal-dirs alone, no
env marker, no role test — per item 68's own 2026-09-22 amendment.

My own verification, both directions, because a hook that blocks worktree commits
would stop every agent in this project:

    # cwd = /home/utig5/dliu/EQdyna (the main checkout)
    pre-commit: REFUSED -- this is the MAIN CHECKOUT, not a worktree.
      --git-dir        = /home/utig5/dliu/EQdyna/.git
      --git-common-dir = /home/utig5/dliu/EQdyna/.git

    # cwd = a linked worktree
    git -c core.hooksPath=testsys/hooks commit --allow-empty  ->  exit 0

It is inert until `install-eqdyna.sh` runs in a given clone (git installs no
hooks on clone and `core.hooksPath` is unset today), so it breaks no session
that is already running.

**Gate.** Both branches were cut from `894cdc1`, which was still `origin/master`
at merge time — no stale base, and the shared-file diffs (`run_e2e.py`,
`run_e2e_full.py`, `install-eqdyna.sh`, `.gitignore`) contain only the intended
additions with no reverted lines. `python3 testsys/run.py unit regression` in the
integration worktree: `SUCCESS unit (exit 0)`, `SUCCESS regression (exit 0)`,
twice — before and after the version bump. The e2e sweep was NOT re-run for this
patch and nothing in the release notes claims one.

**A rule-19 violation, mine to report.** BOTH missions edited `PROJECT_RULES.md`
and `pathway_forward.md` themselves — item 70's rewrote rule 21a's tier paragraph
and its own board row. Those two files have exactly one author and it is
`zofia-kaminska`. I reverted both files to `894cdc1` before landing (integration
commit `02e620d`) so nothing of theirs reached master, and routed the content to
her with the evidence. The generalisable half: my briefs said the rule book was
mine to flip but did not FORBID the edit, and two independent agents did it
anyway on the same night.

## 2. Item 73 — 16 and 32 ranks, the last open question on the jax-MPI A/B

`testsys/perf/run_jaxmpi_ab.py`, `test.tpv104`, five repetitions, arm order
shuffled (base,a / a,base / base,a / a,base / base,a), one cpu set per rank count
per invocation, every ratio paired inside its own repetition,
`--exclude-cpus 0,1,16,17,18,19` (an exclusion by MEASURED effective cores — the
`/proc/stat` probe reads those cpus idle and they then deliver a quarter core).

**Arm `b` was not run, and that is a correction to the board's own unrun
command.** Change b landed at `663d758`, so master's `driver.py` is now IDENTICAL
to the merged branch's — `git diff master mira/jaxmpi-merged-2026-09-22 --
src/python/eqdyna/driver.py` is empty. Running `--arms base,b` today would
measure a file against itself. The only difference left between the arms is
`MPI4NodalQuant.exchange`, which is the non-blocking halo, which is the wave
argument, which is the actual open question.

| ranks | paired speedup, arm a vs base | accepted pairs |
|---|---|---|
| 16 | 1.078, 1.024, 1.077, 0.991, 0.998 — mean **1.034**, sd 0.042 | 5 of 5 |
| 32 | 1.004, 1.026, 1.003 — mean **1.011**, sd 0.013 | 3 of 5 |

The two discarded 32-rank points were discarded correctly: the base arm read min
`EFFECTIVE_CORES` 0.83 and 0.86 while the arm beside it read 1.00, so the ratio
would have compared two different machines.

**Verdict: the wave argument, made at 32 ranks, does not survive above 12 — it
does not appear at 16 either.** The non-blocking halo is now declined on
measurement at 4, 8, 12, 16 and 32 ranks. The robustness argument for it (Irecvs
posted before any Isend cannot deadlock for any neighbour graph, where ascending
`Sendrecv` is deadlock-free only because a slab decomposition happens to give two
neighbours) is untouched by any of this and remains the owner's to weigh.

**A caveat that must travel with these logs.** Absolute per-step numbers on this
box tonight are NOT comparable across repetitions: base-arm per-step at 16 ranks
ranged **84.30 to 110.74 ms** on the same cpu set, and at 32 ranks **32.43 to
57.40 ms**. Only within-repetition paired ratios mean anything, and most of the
variation was our OWN concurrent sessions — one dispatched mission briefly
launched `mpirun -np 16` on this box while testing a mutant that disabled a
refusal, and reported it.

Also worth recording: with the barrier under `prof`, the `exchange` column in the
base arm at 16 ranks reads ~50 ms of a 96 ms step. That is waiting for the
slowest neighbour, not transfer cost — consistent with the known element-count
imbalance, and it is why making the transfers concurrent buys nothing.

Evidence: `docs/NOTES_jaxmpi_1632_2026-09-23.md`, snapshots
`docs/perf_snapshots/jaxmpi_ab_2026-09-2{2_234619,2_235725,3_000153,3_000831,3_001257}.json`,
raw logs `docs/evidence/perf-2026-09-23/wei-perf1632_rep{1..5}.log.gz`, 20 ledger
rows appended (rejected points included, with the rejection on the row).

## 3. Item 64 — setup rewrite, Stage 0 only

New probe, `testsys/perf/run_setup_probe.py`: `eqdyna3d.build_solver_state`
ALONE — no time loop, no MPI, no jax — at 1 process and at 32 concurrent, three
repetitions, `test.tpv104`, effective cores taken from each child's OWN
`getrusage` after the build rather than from a probe before launch, each child
writing its own JSON the moment it finishes (rule 20a).

    n=1   build_s min/mean/max 4.42/4.42/4.42   eff 1.00   maxrss 1650 MB   wave wall  4.6 s
    n=32  build_s min/mean/max 5.73/7.92/9.59   eff 1.00   maxrss 1649 MB each (sum 52.8 GB)   wave wall 9.7 s
    SUMMARY n=1   build_s mean 4.42 over 3 process-samples
    SUMMARY n=32  build_s mean 7.92 over 96 process-samples

**What this settles:**

1. **The "~100 s setup" figure is refuted.** It is 4.42 s serial and 7.92 s mean
   under 32-way concurrency. The figure had never been measured by anything —
   `writeCompTime` is declared `0` at `globalvar.f90:109` and set to 1 nowhere,
   so `compTimeInSeconds(1)` has never printed — and it is now superseded by a
   number from a different instrument rather than merely unsupported.
2. **The MEMORY argument is the one that survives, and it is linear:** 1.65 GB
   per rank, 52.8 GB summed at 32 ranks. On this 1 TB box that is affordable;
   it is what caps rank count on a smaller machine.
3. **The TIME argument is real but small, and only for short runs.** 32-way
   concurrency inflates per-process setup 1.79x (4.42 -> 7.92 s) and the whole
   wave completes in 9.7 s. Item 64's own gap — 32-rank WALL per-step 121.73 ms
   against a differenced 45.76 ms — implies a fixed term of roughly 12 s at
   n=160, and the measured 9.7 s setup wave accounts for most of it. So the
   suspected cause IS the cause, and it is worth ~10 s per run, not ~100 s.

**Stages 1-3 NOT started**, per instruction. What Stage 0 says about them: a
rank-local build would return ~10 s of fixed cost per run and ~51 GB of peak
memory at 32 ranks, and would not move the differenced per-step number the perf
tier gates on at all.

**A methodological finding that came free, and it is the most transferable thing
here.** `EFFECTIVE_CORES >= 0.99` from rusage does NOT mean comparable per-core
throughput. The same `build_solver_state` on the same case measured **19.03 s at
eff 1.00** pinned to cpu 0 while a 16-rank MPI job measured elsewhere on the box,
and **4.42 s at eff 1.00** on a quiet box — 4.3x apart at the same "effective
cores". rusage measures whether the process got a core, never how fast that core
ran. Every ledger row that cites effective cores as its quality gate inherits
this; routed to `zofia-kaminska` as a proposed rule.

Evidence: `docs/perf_snapshots/setup_probe_2026-09-23_000630.json`,
`docs/evidence/perf-2026-09-23/wei-setup0_stage0.log.gz` and
`wei-setup0_smoke.log.gz` (the 19.03 s contaminated point is in the latter, kept
deliberately — it is the evidence for the caveat above).

## 4. Order of operations, and why

The box was idle at 12/64 when this session started, which is the condition every
perf item had been waiting on. Both heavy jobs were run detached and polled by
artifact (rule 20), and they were CHAINED rather than stacked (item 42): the
Stage 0 probe waited on `ALLDONE` from the A/B before launching its first 32-way
wave. Light code missions ran alongside; heavy measurement never overlapped
heavy measurement. The two 32-rank base points that were rejected are the price
of running code missions beside a measurement, and the filter caught both.

## 5. Open, for the next session

- Stages 1-3 of the setup rewrite: decided but not started, awaiting the owner on
  Stage 0's numbers above.
- Item 73's row wants closing by Zofia with the table in section 2.
- `testsys/perf/run_perf.py:172` and `testsys/perf/run_jaxmpi_ab.py:88` rebuild
  fixed in-repo paths and are the same defect class item 70 just fixed;
  `runlock.acquire` is reusable for both. Not queued.
- Tag `v5.16.1` is gated on CI being green (rule 15a) — placed later in the
  session, see section 8.

## 6. Second landing wave — items 74, 75, 77, 80 and two CI guards

All in `v5.16.1` (the tag was held until the whole wave was in), each gated by me
with `python3 testsys/run.py unit regression` exit 0 in an integration worktree,
each branch cut from the then-current master.

| item | what landed | guard |
|---|---|---|
| 74 | `runlock` takes a relative PATH; `run_perf.py` locks `testsys/perf/perf_case` inside `build_perf_case()`, `run_jaxmpi_ab.py` locks `src/python/eqdyna` | `test_perf_tool_locks.py`, 10 checks |
| 75 | `pre-commit` also refuses a commit staging `PROJECT_RULES.md`/`pathway_forward.md` together with any other file | `test_precommit_board_separation_guard.py` |
| 77 | locks for `run_scaling.build_py_case` and `run_numa_scaling.build_case`; `run_jaxmpi_ab` REFUSES when `$EQDYNAROOT` names another checkout | same guard, 10 → 17 checks |
| 80 | `check_board_separation.py` wired as the first step of `test.yml`'s `unit-regression` job | `test_ci_board_separation_step.py` |
| — | `publish.yml`'s `fetch-depth: 0` (see section 7) | `test_publish_image_fetch_depth.py` |

Two things worth keeping out of the table:

**`run_jaxmpi_ab.py` locks the PACKAGE, not the `__pycache__` the board row named.**
Its `stage()` copies arm-specific `driver.py` and `MPI4NodalQuant.py` into
`src/python/eqdyna` and every timing is taken against whatever those files then
hold. Two concurrent invocations in one checkout do not crash — they produce **a
number recorded under the wrong arm label**, which the tool's own docstring calls
worse than no number. I ran five repetitions of that tool tonight; this is the
collision that would have quietly corrupted section 2.

**My own verification of the board-separation hook, run fresh:** staging
`pathway_forward.md` together with `README.md` → `pre-commit: REFUSED -- this
commit MIXES the rule book / board with other files`, exit 1; the same commit
board-only → exit 0. And the measured limits, which matter more than the feature:
git does NOT run `pre-commit` for an automatic merge commit (so the merge path is
unaffected and must stay so), and `git commit --amend` is blind to the check
because `git diff --cached` compares against the commit being replaced.
`check_board_separation.py` closes the amend hole over finished history, which is
what item 80 wired into CI.

## 7. The Docker publish workflow was red on the v5.16.0 tag, and nobody read it

Raised by the coordinator, re-derived by me from the failing run rather than
accepted: run 35819232914 on `894cdc1` — the exact SHA `v5.16.0` is tagged at —
failed two regression guards INSIDE the built image, `test_history_table.py`
("expected 8 bulk-imported commits stamped 2018-11-26, found 0") and
`test_pretag_ci_negative.py`, whose own message named the fix: "shallow or
incomplete checkout: 5 referenced commit(s) do not resolve ... set `fetch-depth:
0`". Cause: `publish.yml` used `actions/checkout@v2` at the default depth 1 and
`Dockerfile:33` is `COPY . /opt/eqdyna`, so the SHALLOW `.git` travels into the
image and both guards read history that is not there.

**The guards were right and the environment was wrong**, so the fix is
`fetch-depth: 0` (`a99902f`) and NOT teaching the guards to skip — a guard that
skips when its evidence is missing is not a gate.

**Consequence, from the run's own step list: `Push image -> skipped`,
`verify-published-image -> skipped`.** The gate runs before the push, so
**`ghcr.io/eqdyna/eqdyna:v5.16.0` was never published and `:latest` was never
moved to it.** Nothing in this project reads that workflow, so a release shipped
without its image and every other gate stayed green. (Registry not queried
directly: this token lacks `read:packages`. The evidence is the job, not the
registry.)

**Rule 15a's pre-tag gate reads ONE workflow's conclusion.** v5.16.0 passed it
honestly — both "Automatic Testing of EQdyna" runs on `894cdc1` succeeded — while
a second workflow on the same SHA was red. Whether 15a should require every
workflow for the SHA is a methodology change at the release boundary; recorded
for the owner, not decided.

## 8. The tag, and proving the fix before the release rather than with it

The coordinator's catch, and it was the right call: `publish.yml` fires only on a
tag push or a dispatch, so tagging would have made the release the experiment.
The workflow's own header says `workflow_dispatch` runs build+gate WITHOUT
publishing. So, in order:

1. `gh workflow run publish.yml --ref master` on `4ee171b` — run 35828671705,
   **success**, nothing published. The shallow-clone fix is proven.
2. `python3 testsys/regression/check_pretag_ci.py --pre-tag 4ee171b` — `PASS  a
   completed, successful Automatic Testing of EQdyna run exists`, exit 0.
3. `git tag -a v5.16.1 4ee171b` and push.
4. The tag's publish run 35828977752: `Gate -> success`, **`Push image ->
   success`**, and the separate `verify-published-image` job pulled the image
   fresh on a clean runner and re-ran the whole gate against it — **success**.

`:latest` now points at v5.16.1. **v5.16.0 still has no image**; backfilling it
is the owner's call, not mine.

## 9. One error of my own

I dispatched a general-purpose agent with the literal prompt "placeholder",
having reached for a continue-an-existing-agent facility this deployment does not
give me. It did nothing (0 tool uses) and the main checkout was verified clean
afterwards, but a full-tool agent with an ambiguous prompt and no worktree
isolation is exactly the containment failure I am here to prevent. Cost: one
wasted dispatch. The rule it implies is that a dispatch is made only with its
isolation and its scope already written.

## 10. Still open after this session

- **Stages 1-3 of the setup rewrite (item 64)** — with the owner, on section 3's
  numbers. Not started.
- Item 64's residual question: nothing has ruled out process-launch/import
  overhead as part of the 32-rank fixed term; the probe measures
  `build_solver_state` only (import is 0.1 s of it).
- `testsys/parity/probe_plastic_traction.py:80` rebuilds a fixed in-repo path and
  is the last unlocked member of item 70's class.
- Four unguarded settings in `.github/workflows/` (publish.yml's `packages:
  write`, the in-image gate's contents, the two `if: github.event_name == 'push'`
  conditions, and the `@v4`/`@v2` checkout divergence).
- No `.dockerignore`: `COPY . /opt/eqdyna` copies whatever is in the build
  context. **A naive fix that excludes `.git` would re-break the two in-image
  history guards** — the trap is worth writing down before anyone takes that row.
- Whether rule 15a should read every workflow for the SHA (section 7).

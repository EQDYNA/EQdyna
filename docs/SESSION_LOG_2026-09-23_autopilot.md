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

---

# Continuation — 2026-09-23 02:20 onward (wei-lin, second conductor)

Resumed at `origin/master` = `5c91598`, clean tree, 0 unpushed, VERSION 5.16.1,
box at 12/64. ~6 h remaining of a 12 h budget granted 20:08, under a standing
ROLLING RENEWAL order: pace for continuity, not for a deadline.

**Grant, stated back.** Minor and patch tags on `master` of THIS repo; merge and
tag authority on `master` directly. No major bump. No force-update or rewrite of
an existing tag. No publish to a package index. Nothing outward-facing beyond
this repo. That last clause is the one that binds section C below.

## A. Item 83 — CLOSED, on commands I ran myself

The row said `v5.16.1` was tagged with no GitHub Release object and that the
regression tier was RED at master HEAD because of it. Both halves are now false.

    $ gh release list --limit 4
    v5.16.1 -- concurrency enforcement   Latest   v5.16.1   2026-09-23T07:16:39Z
    v5.16.0                                       v5.16.0   2026-09-23T04:39:46Z

    $ python3 testsys/regression/test_release_complete.py      # at 5c91598
      PASS  a completed, successful Automatic Testing of EQdyna run exists for
            tagged sha 4ee171b4bc70abcf05f9d524af5fcfb903f167c2 (v5.16.1)
      PASS  v5.16.1 is pushed to origin
      PASS  GitHub Release exists for v5.16.1
    SUCCESS test_release_complete

    $ python3 testsys/run.py unit regression                   # at 5c91598
    SUCCESS unit (exit 0)
    SUCCESS regression (exit 0)                                # 79.8 s wall

Routed to `zofia-kaminska` with that literal output for the row write, together
with the rule the incident paid for: the Release object and the tag push are ONE
action, because a regression guard reads the Release. The cost of separating
them was ~20 minutes of a red tier at master HEAD for a reason unrelated to any
code under test, which is how a real red later gets discounted.

## B. One claim from the previous session, softened rather than repeated

"Re-verified from the registry" is NOT reproducible here and does not stand. The
package is private, anonymous manifest GETs on `v5.16.1`, `v5.16.0` and `latest`
all return 403, and the available token lacks `read:packages`. What the evidence
actually proves is a WORKFLOW conclusion: publish run 35828977752's `Push image`
step succeeded, and the separate `verify-published-image` job pulled the image on
a clean runner and re-ran the whole in-image gate against it — both green. That
is strong, and it is not a registry query. To make the registry claim
reproducible later, a token with `read:packages` is what it would take; recorded
rather than asserted.

## C. The v5.16.0 image backfill — STOPPED, and it is not a re-run

Routed to me as ordinary work on the reasoning that it is a re-run of a
now-fixed workflow. **The code says otherwise, unambiguously, and the
instruction I was given names this exact stop condition.** Three independent
blockers, each sufficient on its own:

1. **A dispatch cannot publish anything.** `publish.yml:42-46` sets
   `TAG=dispatch-${GITHUB_SHA::12}` on any non-`push` event, and `:77-81` gate
   the `Push image` step on `github.event_name == 'push'`. A
   `workflow_dispatch` on the `v5.16.0` ref builds and gates and pushes
   nothing — by design, and I am relying on that same design to gate item 79.
2. **A tag push cannot publish `v5.16.0` without moving `:latest` backwards.**
   `:79-81` push `${TAG}` and `:latest` unconditionally in the same step, and
   `:51-53` re-tag `:latest` on every push event. There is no code path that
   publishes a version tag alone. The owner's instruction was explicit that
   `:latest` belongs to v5.16.1 and must not move back.
3. **Reaching the publish path at all would mean re-pushing an existing tag**,
   which my grant excludes by name — and the workflow file *at that ref* is the
   broken one anyway: `git show v5.16.0:.github/workflows/publish.yml` has
   `actions/checkout@v2` with no `fetch-depth`, because the fix `a99902f` is
   not an ancestor of `v5.16.0` (`git merge-base --is-ancestor` → NO). A
   dispatch on that ref runs the RED workflow, not the fixed one.

So backfilling `v5.16.0` requires a `publish.yml` change — a dispatch input that
publishes a caller-named version and skips `:latest`. That change deliberately
inverts the premise item 82(c) exists to guard ("invert either and a manual
dispatch publishes"), so it is a release-boundary methodology decision, not
ordinary work. **Not attempted. Owner's call.** The cheaper alternative worth
putting beside it: leave `v5.16.0` imageless and record it, since the release
itself shipped and every other gate was honest.

## D. Worktrees reaped: 12 → 7, after checking each for unlanded work

Nine were clean or held only duplicates of evidence already in `origin/master`,
which I verified file by file rather than by eyeballing status — every artifact
the previous session's sections 2 and 3 cite (`testsys/perf/run_setup_probe.py`,
`docs/NOTES_jaxmpi_1632_2026-09-23.md`, the five `jaxmpi_ab_2026-09-2*.json`
snapshots, the seven `docs/evidence/perf-2026-09-23/*.log.gz`, and the ledger at
275 rows against the worktree's own 275) is PRESENT in master. The dirty content
left in `wei-perf1632`/`wei-setup0` was raw `rep*.log` duplicated by the
committed `.gz`, per-child probe JSON superseded by the committed aggregate, and
`nohup`/launcher scratch. Called scratch deliberately (rule 8), not swept.

**Two kept, both owner-held and both locked:** `item32-dx250-refinement`
(`23971a8`, item 32) and `mira/jaxmpi-merged-2026-09-22` (`99a3f18`, ahead 3 of
master — and pushed, `git ls-remote` confirms the ref at that exact SHA, so
nothing there depends on the worktree surviving).

## E. Master HEAD has no CI run, and that is CORRECT — do not chase it

`5c91598` has no "Automatic Testing" run and one is not coming. The four commits
since `d6d846f` touch only `PROJECT_RULES.md`, `pathway_forward.md`,
`docs/NOTES_item70_lock_2026-09-23.md` and `docs/SESSION_LOG_...md` — every one
of them inside `test.yml`'s `paths-ignore`. Writing it down because the absence
looks exactly like a missed trigger, and because it has a real consequence:
**rule 15a's guard requires a successful run for the TAGGED sha, so a tag must
never be cut on a docs-only HEAD.** Every mission in flight below touches
non-ignored files, so this does not bind today.

## F. Landings — six commits, no reverts

| commit | item | what | my own gate |
|---|---|---|---|
| `32fba21` | 83 | Release object confirmed, tier green; rules 15c, 3a | `test_release_complete.py` 8/8 PASS |
| `07b0dca` | 79, 82(d) | root `.dockerignore`; `checkout@v2` → `@v4` | mutation: `.git*` takes the guard RED naming every probe path |
| `73ac3a5` | 81 | `runlock` on `probe_plastic_traction.build_case()` | `SUCCESS test_perf_tool_locks (18 checks)` |
| `d9c9611` | 67 | station header states 8 or 11 per friclaw branch | my own e2e cell + `8 / 8 / 8` |
| `650bd1c` | 84 | `persist-credentials: false` + in-image credential gate | my own mutation, exit 1 / exit 0 |
| `20d8f23` | — | board closures, rule 14a | tier green |

Every branch came back based on a commit master had already moved past, and every
one of them would have REVERTED something if copied wholesale — Zofia's would
have dropped 111 lines of this log, three others would have undone
`.dockerignore` or item 67's Fortran change. None did, because the base diff is
checked before anything is copied. Worth stating as a rate rather than an
anecdote: **6 of 6 returning branches were stale.** On a tree moving this fast
that is the normal case, not the exception.

## G. Item 84 — a GITHUB_TOKEN was being baked into every published image

Found by `haruto-nakamura` incidentally while bumping `checkout@v2` → `@v4`, and
it is the most consequential thing this session produced. Every link verified by
me before I acted on it, because a security claim from a subagent is still a
hypothesis:

1. `actions/checkout@v4`'s `action.yml` declares `persist-credentials: default:
   true` — I fetched the file from `raw.githubusercontent.com` and read it.
2. That default writes the live `GITHUB_TOKEN` into the workspace `.git/config`
   as `http.<origin>/.extraheader`.
3. The same `action.yml` declares `post: dist/index.js`. Cleanup is a POST step,
   so it runs AFTER every main step — including `docker build`.
4. `Dockerfile:33` is `COPY . /opt/eqdyna`, and `.dockerignore` deliberately
   keeps `.git` so the two in-image history guards can read it.

So the token was in `.git/config` at build time and is in a layer of every image
this workflow has pushed, `v5.16.1` included. **Short-lived** — it expires when
the job ends — so this is a credential-in-artifact leak, not a live-token leak,
and the mechanism is unchanged since checkout v2.0.0, so the `@v4` bump did not
introduce it.

Fixed in two layers, because the setting alone is item 82's "green until a
release" shape all over again: `persist-credentials: false` pinned structurally
by `test_publish_image_fetch_depth.py`, and `testsys/check_git_config_no_credential.py`
run as the FIRST line of both in-image gate blocks, before any pip install, so a
regression fails the job instead of publishing quietly. My own mutation of it:

    FAIL check_git_config_no_credential: 2 credential marker(s) found ... (43 line(s) scanned)
      line 43: matched 'extraheader' -- extraheader = AUTHORIZATION: basic ...
    PASS check_git_config_no_credential: ... scanned (42 line(s)), no credential marker among 10 checked

Both verdicts print the line count they scanned. That is the papercuts rule — a
green that never opened the file is this project's recurring failure — and it is
why the refusal path exits 2 on a missing file rather than reporting clean.

End-to-end proof, not just a static one: dispatch run **35834959371** on
`650bd1c` built the image and ran the gate with that check first —
`Gate the image ... -> success`, and `Push image -> skipped`.

**OWNER'S, not mine:** whether the already-published images (v5.16.1 and
earlier) should be deleted and re-pushed. That is outward-facing and outside the
grant. Board row 84 carries it.

## H. Two board rows had an evidence command that could never go green

Item 71's command greps a file the fix was deliberately NOT going to touch — the
coverage landed as a sibling guard (`test_fractal_fault_geometry_derivatives.py`,
`2590c65`, 12 hits, 6 checks green, mutation proofs baked in), so the row's
command reads 0 before the fix and 0 after. It cost a full dispatch, which the
agent correctly spent refusing to duplicate work.

Item 67's command says to find the station file by its `# location = on fault`
header. No file has that line: `stLocStamp` is computed at
`src/fortran/library_output.f90:55` and never written. The command was
unrunnable from birth.

Both rows carried a recent `last checked` date. A date beside a command that
cannot distinguish done from not-done is worse than a blank one, because it
reads as a board being maintained. Now rule **14a**, and pathway row **85**
carries the `stLocStamp` dead variable — same shape as `writeCompTime`, which
this project already knows hid a 511x error.

## I. Nine landings, and the stale-base rate is the finding

`32fba21` (83) · `07b0dca` (79, 82d) · `73ac3a5` (81) · `d9c9611` (67) ·
`650bd1c` (**84**) · `20d8f23`, `ef7b444`, `7ea2f27` (board/rules) · `4c9cc26`
(CLAUDE.md) · `bef0f05` (82abc). No reverts.

**Every single returning branch — nine of nine — was based on a commit master
had already moved past**, and every one would have reverted something had it
been merged wholesale: Zofia's first pass would have dropped 111 lines of this
log, Iris's item-79 branch would have undone two rules, the credential branch
would have undone item 67's Fortran change and the `.dockerignore`, the item-82
branch would have undone the CLAUDE.md fix. None did, because gate axis 4 runs
before anything is copied: `git diff --name-only origin/master <branch>`, take
only the intended files, then re-read the removed lines.

On a tree moving at nine commits in four hours this is the NORMAL case. Stop
treating a stale base as an agent error — it is a property of dispatching into
worktrees at all, and the only defence is the diff, every time, with no
exceptions for "it's just a test file".

One agent caught its own version of this mid-mission: it diffed against the
local `master` ref, which still pointed at `5c91598`, noticed the numbers were
impossible, and re-ran everything against `origin/master`. **Local `master` in
a worktree is not `origin/master`** and a `git fetch` does not move it.

## J. Item 82 is closed on all four sub-items, and (a) is the transferable one

`testsys/regression/test_publish_image_fetch_depth.py` went 5 → 13 printed
checks. Two of them are worth naming because a reasonable guard would have
missed both:

- **The two in-image gate blocks are pinned as EQUIVALENT to each other**, not
  merely each present. The risk was never a missing block; it was an edit to one
  and not the other. Mutation: remove `pip3 install jax` from the
  `verify-published-image` block only → the guard emits a line-by-line drift
  diff. A per-block check passes that mutation.
- **The `:latest` gate has no YAML `if:` at all.** It is a shell `if` inside the
  `Build image` script, so it has to be parsed out of the script body. A guard
  that walked step conditions would have silently skipped the one gate that
  keeps `:latest` from moving backwards — which is the property I relied on
  twice tonight when I fired dispatch runs that deliberately published nothing.

And (a), `packages: write`: dropping it changes NO observable behaviour on any
dispatch run. It only breaks a real release. Nothing inferential finds that; it
needs a structural read of `job['permissions']`. That is item 78's shape for the
fourth time tonight — **a premise that is only ever exercised by a release is
invisible to every gate that is not one.**

My own mutation, fresh: flipping `Push image`'s condition to
`workflow_dispatch` →

    FAIL: ... expected `if: github.event_name == 'push'`. Anything else and a
    workflow_dispatch smoke run publishes to ghcr.io.

restore byte-identical, sha256 `02189fb7...cddc5569b`, tier green.

## K. A cost I caused, and the fix belongs on the board

Nine commits in four hours queued nine full CI matrix runs, and **two of them
were for agent BRANCHES whose content was already merged to master**
(`0951165`, `1379647`). `test.yml` triggers on a push to ANY branch, so every
mission that pushes its work — which they must, it is how I fetch it — burns a
seven-job matrix on a commit nobody will ever tag. I cancelled both to free
runner slots for the release gate, which was by then six deep and blocking the
tag.

The mitigation is a branch filter on `test.yml`'s push trigger (master plus
`pull_request`), and it is NOT something to change during a release cut, so it
is recorded rather than done. The general form is worth more than the fix: **a
CI trigger that fires on work nobody will merge converts throughput into
queue**, and it does it invisibly, because every individual run looks correct.

## L. v5.16.2 released — `2964325`, and the release gate earned its keep

| | |
|---|---|
| tag | `v5.16.2` (annotated) at **`2964325`** |
| Release | "v5.16.2 -- credential-in-image fix and publish-path guards", Latest, `2026-09-23T09:54:37Z` |
| rule 15a | run **35839792140**, all 7 jobs success; `check_pretag_ci.py --pre-tag 2964325` → PASS, exit 0 |
| rule 15c | tag push and `gh release create` were one action |
| after | `test_release_complete.py` **8/8 PASS** |
| image | publish run **35845693679** success — `Push image -> success`, and `verify-published-image` pulled the tag fresh on a clean runner and re-ran the whole gate: success |

**The registry question, finally answerable in a reproducible form.** The
package is private and this token lacks `read:packages`, so an anonymous
manifest GET still returns 403 — that has not changed and no log should claim
otherwise. But `verify-published-image` IS a registry round-trip: it pulls the
pushed tag on a fresh runner and re-runs the in-image gate against what comes
back, which for v5.16.2 includes `check_git_config_no_credential.py`. So the
claim that stands is **"the published v5.16.2 image was pulled from the
registry and passed the full gate"**, evidenced by that job. Not "the registry
contents were verified".

## M. The release commit took master red, and that is the gate working

`7144a98` bumped `VERSION` to 5.16.2 and left the runtime banner at
`src/fortran/eqdyna3d.f90:17` saying 5.16.1. Rule 11 already requires those to
move together. `test_version_banner.py` went red, so **the regression tier was
red at master HEAD** — the precise stop-everything condition rule 3a was written
for, seven hours earlier, in this same session. CI run 35838379925 failed on
`unit-regression` and the tag was correctly never cut.

I found it independently, by running the tier myself at master HEAD before
gating the tag rather than trusting that a release commit is safe. Fixed in
`2964325`: banner bumped, rebuilt, tier green, and the `test.tpv8 x fortran`
oracle re-run because the change touched `src/fortran` —
`max|diff|=3.051760e-11` against the `1.0e-08` bound, unchanged.
`grep -rn "Welcome to EQdyna" src/ scripts/` returns exactly one line, so there
was no second copy hiding.

**Two failures, not one, and they need two fixes.**

1. The release commit was pushed without the tier being run on it. The guard
   that catches this already existed and already works; nothing made anyone run
   it. That is procedural and mechanically fixable.
2. **My own brief forbade the fix.** I wrote "do NOT edit anything under
   `src/`" into the release brief, to stop a release turning into a refactor.
   That made rule 11 unsatisfiable — a VERSION bump REQUIRES a `src/` edit. The
   engineer halted and reported rather than guessing which rule to break, which
   was correct, and it cost about forty minutes.

The generalisation is the part worth keeping: **a scope restriction is itself a
rule and can conflict with another rule.** When it does, halting is right
behaviour and the defect is in the restriction. I wrote the restriction; the
cost is mine, not the engineer's.

## N. Also mine: the release agent ended its turn on a wait

It pushed the release commit, reported "CI is in progress, I'm watching it in
the background and will resume once it completes", and exited. It was not
watching anything — it had stopped, and this deployment gives me no facility to
resume an agent. Had I taken that report at face value the release would have
sat parked indefinitely. I drove the tag myself instead.

This is the same failure the previous session recorded against itself and the
same one my own instructions name first. It is evidently not enough to hold
that rule myself: **a dispatched agent whose task ends in a wait will park it**,
so either the brief must forbid ending on a wait, or the waiting must stay with
the conductor. Tonight the right split was the second — the polling was mine to
do all along.

## O. Full local sweep on the released tree — 31/31

`python3 testsys/run.py all` on `e5011fa` (the v5.16.2 tree plus docs-only
commits), box quiet at load 7.7 rising to a peak of 68 during the python cells:

    ran       : 31 of 40 cells (31 passed, 0 failed)
    not gated : 9 declared-unsupported (every case x python-jax-mpi except test.tpv8)
    wall clock: 1321.6s
    SUCCESS unit / SUCCESS regression / SUCCESS e2e

Evidence `docs/evidence/perf-2026-09-23/wei-s2_sweep_all_v5.16.2.log.gz`,
snapshot `docs/perf_snapshots/e2e_cells_2026-09-23_052913_3839937.json`, 31
ledger rows.

**No minor bump cut, deliberately, against the letter of the version scheme.**
The scheme says minor "after a clean fast/full sweep with accumulated patches".
The patches are already in v5.16.2, tagged ninety minutes earlier; a v5.17.0
whose diff is one evidence file would be a tag nobody can read. The sweep is
VALIDATION of v5.16.2, recorded as such. Stating the deviation rather than
quietly skipping it.

## P. A red master CI run I walked past for two hours

Run **35833383579 on `73ac3a5`** (item 81's landing) **FAILED** on
`unit-regression`, and I did not look at it. I was polling for the runs I cared
about — the tag's — and a red on master slid by underneath, which is exactly
what rule 3a says must stop everything. It did not stop anything, because I
never read it.

**Diagnosis, from the run's own log:** `test_stop_exit_status.py`. The exit CODE
was right — `mpirun -np 2` on an empty case exited 21
(`ERR_INPUT_FILE_MISSING`), no hang — and the registry (22 codes, unique,
1-125) and the README exit-code table both matched. What failed was
`the refused run printed no FATAL block`: the FATAL text was not captured. A
two-rank MPI stdout capture race on a shared runner, not a solver defect.

**It did not recur:** `d9c9611`, `4c9cc26`, `bef0f05` and `2964325` all ran the
same job green afterwards, and the tagged SHA's run is 7/7. So v5.16.2 is not
in question. But **an intermittent gate is not a gate** — a guard that fails one
run in five teaches readers to re-run instead of read, which is the same
discount rule 3a exists to prevent. Routed as a board row.

Two lessons, and the first is mine: **poll the run list, not the run you want.**
Filtering `gh run list` down to the SHA I was gating hid a red on the same
branch. And the reason it was survivable is worth naming — I ran
`unit regression` locally at every single landing, so the local gate had already
cleared 73ac3a5; the CI red was the flake, not the code.

## Q. Item 87's premise, measured — and my own framing in section K was WRONG

I wrote that CI firing on branch pushes is "invisible waste" and proposed a
`branches: [master]` filter. Measured over the last 200 runs
(`gh run list --workflow=test.yml`), push-triggered only:

| | runs | success | failure | cancelled |
|---|---|---|---|---|
| master | 146 | 129 | 16 | 1 |
| non-master | 46 | 35 | **7** | 3 |

Seven non-master reds is not zero catch-rate. And the identity of those seven is
decisive: **six of them are TAG pushes** — `v5.16.1`, `v5.8.2`, `v5.7.0`,
`v5.6.1`, `v5.6.0`, `v5.3.6` (a tag push reports the tag as `headBranch`). Only
ONE was a real feature branch (`iris/perf-ledger-platform-and-e2e-capture`).

So `on: push: branches: [master]` would **stop testing tag pushes entirely**,
because naming `branches:` without `tags:` drops tag events. That is a gate
removal, not waste trimming, and it would remove the signal that caught six
reds. My proposed fix was wrong; the row must carry the measurement and not the
proposal.

This is the house rule about unmeasured signal removal, applied to me: I wrote
a removal proposal into a board row on a plausible mechanism and the count
refuted it inside ten minutes. Falsify against the outcome, not against the
defect you predicted.

**And a second finding falls out of the same table:** `v5.16.1`'s tag-triggered
run on `4ee171b` was RED, and nothing read it — `check_pretag_ci.py`'s
`drop_tag_triggered_runs` excludes tag-triggered runs by design (that exclusion
IS rule 15a's mechanism), so the pre-tag gate was honestly green while the
tag's own run failed. Near-certainly the ordering incident of section M's
ancestor: the run started before `gh release create`, so
`test_release_complete` saw no Release. **v5.16.2 does not repeat it** — its tag
run 35845693484 has `unit-regression: success`, so the Release existed by the
time the guard ran. Rule 15c worked. But nobody was reading tag-run conclusions
at all, and that is item 78's question ("should 15a read every workflow for the
SHA") arriving from a third direction.

## R. Item 84 confirmed from our OWN CI log, not just from the action's source

The post-job cleanup of run 35833383579 prints it literally:

    [command]/usr/bin/git config --local --name-only --get-regexp http\.https\:\/\/github\.com\/\.extraheader
    http.https://github.com/.extraheader
    [command]/usr/bin/git config --local --unset-all http.https://github.com/.extraheader

The credential WAS in `.git/config`, and the unset happens in **Post job
cleanup** — after every main step. That is item 84's whole mechanism, observed
in this repository's own run rather than inferred from `action.yml`. Worth
having: the fix shipped on a chain of documentation reads, and this closes it
with direct evidence.

## S. Item 85 — the discarded value was also WRONG, which is the point

Landed `db0b006`. `stLocStamp` was not merely computed and dropped; the value it
computed was incorrect. The header built its down-dip field from
`abs(xonfs(2,...)/1000.d0)` — raw VERTICAL depth — and labelled it
`km down-dip`. **Six lines above, `library_output.f90:48`, the FILENAME field
already divided by `dsin(fltxyz(2,4,1))`.** So on any dipping fault the filename
and the header of the same file disagreed about the same station.

Verified myself from three independent points rather than from the report:
`meshgen.f90:905` matches `xonfs(2,...)` against `nodeCoor(3)` (vertical z);
`meshgen.f90:917` is `un(3,...) = dcos(fltxyz(2,4,*))`, so that slot is the dip;
and `:48` already encodes depth/sin(dip).

**The mission's own e2e evidence proved nothing about its own change, and this
is the merge-gate axis that catches it.** It ran `test.tpv8`, which is dip 90,
where `sin = 1` and the fix is arithmetically a no-op — a case that falls
through the new path rather than triggering it. I ran `test.tpv36` instead
(`tpv36_37_common.py:16`, `par.dip = 15`):

    test.tpv36  fortran  SUCCESS  174.9s  max|diff|=0.000000e+00 bound=1.0e-06, 3477 nodes
    test/test.tpv36/faultst000dp010.txt:  # location = on fault, 0.0 km along strike, 1.0 km down-dip
    test/test.tpv36/faultst000dp030.txt:  # location = on fault, 0.0 km along strike, 3.0 km down-dip

`bStations.txt` gives that first station a depth of `-0.25881904510252074` km,
which is exactly `1.0 x sin(15 deg)`. So the header used to read **0.3 km** in a
file named `dp010` and now reads **1.0 km**, agreeing with its own filename for
the first time.

**The lesson is sharper than the row was:** a value that is computed and
discarded is not dormant, it is UNVERIFIED. This one had been wrong for as long
as it had been unread, and the discard is precisely what hid it — the same
shape as `writeCompTime`, which this project already records as having hidden a
511x error. Treat every write-only variable as a defect that has not been
looked at yet.

## T. A third process failure, and it is mine again

The final board pass was pushed **straight to master without passing through my
gate**. It was fine — two board files, tier green at `7587244`, verified
retroactively — but it landed unreviewed.

The cause is my wording. Earlier briefs ended "commit to your branch and push
it; do not merge." The last one ended "Push; do not merge, do not tag", and
"push" with no object reasonably reads as push to master, with "do not merge"
covering only branch merges. **Same defect class as the `src/` restriction that
halted the release engineer** (section M): I wrote the constraint, the agent
read it correctly, and the constraint was wrong. Three of tonight's four
process failures originate in brief wording, not in agent behaviour.

The fix is mechanical and belongs in the brief template rather than in anyone's
memory: state the ref explicitly — "push to `refs/heads/<your-branch>`; do NOT
push to master" — every time.

## U. Close-out

| | |
|---|---|
| tag | `v5.16.2` at `2964325`; Release Latest; tag-triggered CI run **35845693484 success**; publish run 35845693679 success incl. the fresh-pull re-gate |
| master | `7587244` at close, clean, level with origin |
| landings | 13 commits, **0 reverts** |
| sweep | 31/31 cells, all three tiers, on the released tree |
| worktrees | 4 → primary, `wei-s2` (mine), and 2 owner-held LOCKED (`item32-dx250-refinement`, `mira/jaxmpi-merged-2026-09-22`, the latter pushed to origin at `99a3f18`) |
| CI at close | `db0b006` 6/7 green (`unit-regression` and both fortran e2e jobs included); only `e2e-ci-python-cheap` outstanding, and this landing touches no Python |

The primary checkout was fast-forwarded from `5c91598` to current master.
It had been 12 commits behind all session, and a stale local `master` already
misled one agent tonight into diffing against the wrong base — leaving it
behind would have handed that same trap to the next session.

---

# Continuation — 2026-09-23 06:25 onward (wei-lin, third conductor)

Resumed at `origin/master` = `078bb4f`, clean tree, 0 unpushed, VERSION 5.16.2,
tag `v5.16.2` at `2964325`, 4 worktrees, box at load 7.5.

**Grant, stated back.** Minor and patch tags on `master` of THIS repo; merge and
tag authority on `master` directly. No major bump. No force-update or rewrite of
an existing tag. No publish to a package index. Nothing outward-facing beyond
this repo. Item 86 (expired token in already-published images), the v5.16.0
image backfill, item 64 Stages 1-3, and items 63, 56, 57, 17, 19(b), 32 and rule
15b's staleness bound are OWNER-HELD and not decided here.

## AA. Item 33 — the jax-vs-Fortran ms/step table, on the quietest box this item has ever had

The precondition this row has carried since 2026-09-16 ("REQUIRES AN IDLE BOX")
and that the 2026-09-21 pass declared UNSATISFIABLE (0 of 64 cpus under the
strict 0.2 ceiling, twice) was **satisfiable this morning and I measured it
before doing anything else**: at 06:28, **53 of 64 cpus were under 10% busy** and
exactly 7 were pegged at 100% (seven foreign single-core jobs, which is the
entire load average of 7.3). Nothing of mine was running: no agent was dispatched
until all three sweeps were done, because the largest source of variance in the
last two sessions was our own concurrent work.

`testsys/perf/run_scaling.py`, `test.tpv104`, compact placement, per-step by
difference (`n_lo=20`, `n_hi=60`), `--repeats 2` (min within a sweep),
**strict `--busy-ceiling 0.2`, never overridden**, capped at 16 cores, jax only,
no GPU. **Three full independent sweeps**, 06:29-07:07.

| cores | Fortran ms/step (3 sweeps) | jax ms/step | jax / Fortran |
|---|---|---|---|
| 1 | 925.28, 909.05, 916.37 | 607.47, 607.11 | **0.67x** |
| 2 | 464.36, 466.14, 463.58 | 405.47, 405.11, 407.95 | 0.87x |
| 4 | 235.72, 236.52, 238.66 | 358.06, 337.08, 358.01 | 1.43x |
| 8 | 122.94, 116.38, 117.35 | 268.78, 319.16 | 2.31x |
| 16 | 55.97, 56.18, 60.14 | 234.44, 247.88 | **4.19x** |

(ratio of the fastest point at each core count; each cell's own spread is printed
beside it rather than averaged away.) Self-relative scaling from each engine's
own 1-core point: **Fortran 1.00 / 1.96 / 3.86 / 7.81 / 16.24x**, **jax 1.00 /
1.50 / 1.80 / 2.26 / 2.59x**.

Three points are worth more than the table:

1. **jax is FASTER than Fortran on one core — 607 vs 916 ms/step, and again at
   2 cores.** The crossover is between 2 and 4 cores. Every "Fortran is faster"
   statement about this code is a statement about PARALLEL SCALING, not about
   per-core work.
2. **The standing 2026-09-18 Fortran curve was contaminated and is superseded.**
   It recorded 161.69 ms/step at 8 and 119.15 at 16 (speedup 7.95x); today the
   same tool, same case, same method reads 117.35 and 56.2 (**16.24x**, near
   linear). The old numbers were taken at loadavg 24-30. Fortran's scaling was
   never the problem it looked like.
3. **The "4.09x jax at 16 cores" figure does NOT reproduce.** Today jax reaches
   2.59x at 16, which agrees with the OLDER 2.52x-at-16 figure that 4.09x was
   thought to supersede. Both of the higher figures (4.09x, 3.38x) came from the
   same contaminated night, and they were high because the Fortran denominator
   was slow, not because jax was fast. jax's plateau past 8 cores stands exactly
   where the HLO evidence put it (sub-linear `outer_dimension_partitions`,
   1 -> none, 4 -> "2", 8 -> "3", 16 -> "4").

Every point passed the strict per-cpu ceiling on the cpus it actually used; one
configuration was SKIPPED as busy in each sweep (`jax th=1`, then `th=16`, then
`th=8`) and recorded, not overridden — which is why there are three sweeps and
not one: the union covers every cell twice or more. Reproducibility: Fortran
within 3% at 1/2/4, 5.6% at 8, 7.4% at 16; jax within 1% at 1/2, 6% at 4 and 16,
and **19% at 8 cores** (268.78 vs 319.16), which is the one cell I would not
quote to three digits.

Method note carried from the owner: `EFFECTIVE_CORES >= 0.99` is necessary and
NOT sufficient (rule 6a, item 76) — the defence used here is not the rusage
figure, it is three independent sweeps with the per-sweep spread printed and the
box tenancy (7/64) recorded in every snapshot.

GPU was not touched and 32 cores was not measured: both are off this queue by
the owner's instruction.

Evidence, all committed: snapshots
`docs/perf_snapshots/scaling_2026-09-23_0{64050,65349,70641}_item33_jax_vs_fortran.json`,
raw logs `docs/evidence/perf-2026-09-23/wei-s3_item33_rep{1,2,3}.log.gz`,
27 rows appended to `docs/perf_ledger.jsonl` (append-only verified: 27
insertions, 0 deletions against `origin/master`). Routed to `zofia-kaminska` for
item 33's row.

## BB. Three landings, each gated by my own run, 12:14-12:40

| commit | item | what | my own gate, not the mission's |
|---|---|---|---|
| `c4666ae` | 88 | `output_offfault_st`: stamp emitted, `', '` separator added, duplicate legend line dropped | `test.tpv8` writes **11 body station files, two at 12 km depth** — a case that TRIGGERS the path; header now reads `# location = -3.0 km off fault, 12.0 km along strike, 12.0 km depth`, legend 7 lines against `awk NF` = 7; `test.tpv8 x fortran max\|diff\|=3.051760e-11` vs bound `1.0e-08` |
| `32744b2` | 89 | `test_stop_exit_status.py` deflaked: FATAL text asserted where the launcher cannot drop it | 20/20 green; MY OWN mutations — delete the FATAL banner -> RED in probe B (`1342 byte(s), none containing it`), delete `flush(6)` -> RED (`164 byte(s)`); `errorCodes.f90` restored byte-identically, sha256 `f7773da6...b4d85d` both sides |
| `8a3eba8` | refactor round 1 | one copy of the NUMA `{node:[cpu]}` inversion and of the per-step-by-difference arithmetic (4 files, +43/-12) | 22,000 random cases comparing each extracted function against the expression it replaced, repr-exact, 0 differ; a REAL `run_scaling` run exercising both engines, arithmetic re-derived by hand from the printed walls; `SUCCESS test_perf_tool_locks (18 checks)`; unit+regression |

**Item 88's dip question, answered rather than assumed.** No `dsin(dip)`
correction applies off the fault: `readInputFiles.f90:218` reads `x4nds` as a
plain (x,y,z) triple and `meshgen.f90:625` matches its third component against
the raw mesh node `nodeCoor(3)`. A body receiver off the fault plane has no
down-dip coordinate to convert to. I re-derived the seven-column legend myself
from `eqdyna3d.f90:180-190` (station outer, `iDof` 1..3, disp/vel inner) and
`driver.f90:229-231` before accepting the mission's table.

**Item 89's cause is the launcher, and that is worth more than the fix.** Open
MPI discards already-written, already-FLUSHED bytes when `MPI_Abort` tears the
job down: a rank writing 20000 lines returned 2,264,946 of 2,879,970 bytes with
an immediate abort and 2,879,970 of 2,879,970 on all 20 runs with a 1 s sleep
first. Same writer, same capture code, only the abort timing differs. So the
guard now asserts the exit status where the user actually sees it (piped
`mpirun`) and the FATAL text where bytes cannot be lost (each rank redirected to
its own file), and neither claim is softened. Any OTHER guard in `testsys/` that
asserts on text captured through `mpirun` past an `MPI_Abort` has the same
exposure — not swept, recorded.

**Round 1's finding is worth more than its diff.** The guard that protects the
four-times-patched "rmtree without a lock" class (`test_perf_tool_locks.py:644-695`)
asserts the SOURCE SHAPE of three function bodies, and therefore FORBIDS the
deduplication that would prevent a fifth occurrence — proved by mutation, not
argued. A guard written against source text can forbid the fix for the defect it
guards. `iris-vermeulen` is in flight now converting those assertions to
behavioural ones; round 3 does the dedup underneath it, in that order.

**Stale-base rate this session: 3 of 3.** All three branches were cut from
`078bb4f` and would have reverted `b978637` (the item-33 evidence) wholesale.
None did; only the intended files were taken, every time.

**Three merged-branch CI runs cancelled** (`2b4d2d0`, `9229821`, `b5da34a`) with
seven `test.yml` runs in flight at once. Third session running to pay item 87's
cost by hand.

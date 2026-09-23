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

## CC. Landings 4-6, and a brief of mine that conflicted with an agent's own rules

| commit | what | my own gate |
|---|---|---|
| `7484b48` | board pass — items 33, 88, 89 CLOSED; rows 90-95 opened; rules 3c, 4d, 10a | tier green at the merged tree; her branch was stale-based and would have deleted 43 lines of this log, so only the two board files were taken |
| `928c4d4` | item 90 — `test_perf_tool_locks.py` now asserts BEHAVIOUR, not source text | MY OWN mutation: an `rmtree` inserted before the acquire in `run_scaling.build_py_case` takes the guard RED on three checks and prints the sequence `[('rmtree', False), ('rmtree', True), ...]` — destroyed unlocked, then locked. Restored byte-identically (`5a7591fd…9c15`), guard back to 18 checks |
| `11bde96` (branch only) | refactor round 2 — `src/python/eqdyna` clarity, +86/-86 | NOT LANDED YET: held for my own bit-identity re-run behind the measurement window |

**Rule 22 fired for the second time in two sessions, and again the defect was
in my brief.** Round 2 came back with the work complete and UNCOMMITTED: my
brief said "commit and push to `refs/heads/<your-branch>`", and that agent's own
operating instructions forbid it from committing at all. It halted and named
both sides instead of picking one, which is exactly the behaviour I ask for. I
committed it myself from its worktree (`11bde96`, pushed to
`kai/port-clarity-2026-09-23`) and the cost was zero. Worth recording because
the fix is not "tell agents to commit" — it is that a conductor's brief cannot
assume an agent's own constitution permits what the brief asks, and the cheap
form is to say "commit if your own rules allow it; otherwise leave it
uncommitted in the worktree and tell me, and I will".

## DD. The owner re-ordered the queue mid-session, five times, and the table I produced this morning got a correction I did not catch

**The correction, and it is real:** item 33's table compares Fortran's genuine
3D MPI decomposition (`run_scaling.py:169` `DECOMP` — 8:(2,2,2), 16:(4,2,2))
against jax running as ONE PROCESS under `numactl` on N cpus, with XLA and
OpenBLAS auto-sizing their pools from the affinity mask (`run_scaling.py:481`,
`PY_THREADS` against `FORTRAN_RANKS` — the module's own constant names say it).
So "Fortran 16.24x, jax 2.59x" is substantially MPI-vs-threads, not a property
of the port. The ms/step numbers stand; the label was wrong. I verified the
asymmetry in the source myself rather than taking it: `time_one_py` launches a
single `sys.executable -c` under `numactl`, `run_fortran` launches `mpirun -np n`.
The ledger's `ranks` field on a `python-jax` row means CPUS IN ONE PROCESS'S
MASK, which reads as MPI ranks and is actively misleading; the fix is a
`parallelism: threads|mpi` discriminator, not a rename (the ledger is
append-only and a rename breaks every historical row's comparability).

**Queue as of 14:05, owner's order:**
1. 1D jax-MPI table — same shape, jax column from real ranks via
   `testsys/perf/run_mpi_scaling.py` (which already differences each RANK'S OWN
   solve time and REFUSES a non-positive result, `:309-315`). Outside the gated
   matrix; `test.tpv104` is not opted into `matrix.PY_MPI_RANKS` and will not be.
2. Correspondence audit (sophia) + the Fortran-is-gold rule (zofia) — read-only
   and doc-only, dispatched in parallel, explicitly forbidden from running
   anything heavy.
3. Rank-local 3D decomposition (items 43 + 64 Stages 1-3 together), design
   reported BEFORE it is built.
4. Port optimization, starting with the unconfirmed 24.42-vs-11.83 ms/step
   compute inflation.
5. The two remaining refactor rounds, re-scoped from `src/` to `testsys/`.

**The measurement window is a scheduling constraint, not a preference.** The
morning's three sweeps were clean because 53 of 64 cpus were under 10% busy and
NO agent of mine was running. The TPV30 hunt is on the box now (one process at
99.9%, 200 s in), so the 1D table waits for it rather than being measured
through it — the alternative is the 19% spread that already cost sweep 1 its
8-core point.

## EE. Item 33, second table — Fortran 3D MPI vs jax-MPI 1D slab, 1-16 ranks

`testsys/perf/run_mpi_scaling.py --case test.tpv104 --ranks 1,2,4,8,16
--n-lo 20 --n-hi 60 --max-busy 0.2 --syncs halo --warmup --exclude-cpus 0,1`,
two independent sweeps, 08:45-09:05. Per-step by difference over **each RANK'S
OWN SOLVE TIME** (the tool refuses a non-positive result outright,
`run_mpi_scaling.py:309-315`), strict ceiling 0.2 never overridden, cpus 0 and
1 excluded by measurement (they read busy 0.00 and then deliver ~0.39
effective cores). **Run OUTSIDE the gated matrix**: `test.tpv104` is not opted
into `matrix.PY_MPI_RANKS` and this did not change that.

| ranks | Fortran ms/step (2 sweeps) | jax-MPI ms/step | jax/Fortran (best pair) |
|---|---|---|---|
| 1 | 982.19, 924.29 | 763.81, 859.63 | **0.83x** |
| 2 | 472.25, 468.54 | 517.21, 502.50 | 1.07x |
| 4 | 232.81, 235.98 | 249.64, 299.68 | 1.07x |
| 8 | 120.91, 120.06 | 191.23, 187.17 | 1.56x |
| 16 | 59.20, 65.24 | 93.67, 116.69 | **1.58x** |

Self-relative: **Fortran 1.00 / 1.97 / 3.97 / 7.70 / 15.61x**; **jax-MPI 1.00 /
1.52 / 3.06 / 4.08 / 8.16x**.

**Against the morning's threaded table, this is the headline.** At 16 the
threaded jax path read 234-248 ms/step and 4.19x of Fortran; real ranks read
**93.67 ms/step and 1.58x**. Most of what looked like a port deficit was the
measurement comparing a 3D-decomposed MPI solver against threads in one
process.

**Where it turns over, which is what the owner asked for: between 4 and 8
ranks.** The ratio is 1.07x at both 2 and 4 — the port tracks Fortran exactly
— and then goes to 1.56x at 8 and stays there at 16.

**The mechanism is visible in the run's own per-rank output, not inferred.**
Exchange cost per step, 8 ranks: `[54.1, 85.4, 92.7, 51.5, 14.6, 3.5, 9.2,
3.5]` ms against a 191 ms step. At 16 ranks: `[57.2, 57.0, 55.4, 50.9, 48.4,
49.2, 54.5, 56.6, 29.6, 16.7, 50.0, 13.5, 3.0, 3.5, 2.7, 2.8]` against a 94 ms
step. The halo cost per rank does NOT fall as ranks rise — it stays at roughly
50 ms on the interior ranks while the compute half halves — which is exactly
the 1D-slab prediction: the slab's halo SURFACE stops shrinking while its
volume does. That is the case for 3D decomposition, stated in the port's own
numbers rather than from theory.

**Two caveats that travel with this table.**
1. **Reproducibility is worse than the morning's.** Fortran repeats within 1%
   at 2/4/8 and 10% at 16; jax-MPI spreads 12% at 1 rank, 20% at 4 and **25%
   at 16** (93.67 vs 116.69). Both sweeps ran while my own TPV30 hunt held one
   core at 100% — that is a mission of mine, not a foreign tenant, and the
   morning's clean table had none. Every point passed the strict per-cpu
   ceiling on the cpus it used and every rank reported EFFECTIVE_CORES 1.00,
   which per rule 6a is necessary and NOT sufficient. Read the 8- and 16-rank
   ratios as ~1.5-1.6x, not as three digits.
2. **At 1 rank the MPI path is SLOWER than the threaded path on the same
   case** — 763-860 ms/step against the morning's 607 ms/step on one core.
   Same solver, same case; the difference is `driver.run_mpi` plus rank-local
   setup versus `driver.run`. Flagged as an observation to confirm, not a
   finding: the two came from different tools with different placement.

Evidence: `docs/perf_snapshots/mpi_scaling_2026-09-23_08{4522,5541}_tpv104.json`,
`docs/evidence/perf-2026-09-23/wei-s3_jaxmpi1d_rep{1,2}.log.gz`, 20 append-only
ledger rows.

## FF. Landings 7-10, and the brief defect that has now fired twice

| commit | what | my own gate |
|---|---|---|
| `013f762` | rule 23 — Fortran is the reference implementation; the port follows its NUMERICS, not its file layout | counts verified at the merged tree (47 headings, 23 numbered rules); tier green |
| `6c33840` | `docs/fortran_python_correspondence.md`, all 26 `.f90` classified counterpart / folded / absent (rule 23 clause 4) | TWO corrections of my own, below |
| `56d2401` | item 43 + 64 design, DESIGN ONLY | four of its load-bearing claims re-read in the source by me |
| `090d2cf` | refactor round 3 — one shared case-builder (`testsys/perf/perflib.py`) for four tools, net -22 lines | MY OWN mutation on a THIRD tool: delete the acquire in `run_scaling.build_py_case` -> 4 checks RED; restored byte-identically (`2080eeb6…497866`); 18 checks green; real invocation through the shared builder |

**The audit's headline finding was wrong, and checking cost one grep.** It
reported that an inverted element makes Fortran stop while Python produces a
finite wrong answer. Both port determinant paths raise:
`assembleGlobalMass.py:167` (`compute_element_det`) and `:461-464`
(`compute_element_shape`, the path `eqdyna3d.py:301` actually uses). What
survives is narrower and still true — the port has no NUMBERED exit-code
contract. The second correction: `rdampm` appears exactly twice in
`src/fortran` (declaration at `globalvar.f90:179`, use at
`assembleGlobalKU.f90:19`) with no reader, so both sides are 0 for every
runnable case and the "port hardcodes 0.0" gap is a doc fix. A subagent's
report is a hypothesis even when the subagent is auditing someone else's work.

**Round 3's dedup was unblocked by round 1's finding, in the right order.**
Round 1 found that the guard forbade the fix; `iris-vermeulen` made the guard
behavioural (`928c4d4`); round 3 then did the dedup under it. Three missions,
one finding, no step skipped — and the proof that the order mattered is that
the same dedup attempted before `928c4d4` turned the guard RED.

**Rule 22 fired twice today and both times the defect was mine.** Rounds 2 and
3 both came back complete and UNCOMMITTED because my brief said "commit and
push" and that agent's own constitution forbids committing at all. Both halted
and named the conflict instead of picking a side, which is the behaviour I ask
for. The brief template fix: *"commit if your own rules allow it; otherwise
leave the work clean and uncommitted in the worktree and tell me — I will
commit it."* Cost so far: zero, because both said so plainly. It would have
cost the whole mission if either had silently not-committed and reported
success.

## GG. The 3D design, and the question it was dispatched to answer

`56d2401`, design only. **frt canonicalisation stays honest under rank-local
numbering and neither `frt_canonical.py` nor `compare.py` changes** — I
verified the chain myself rather than accepting it: `write_frt` builds rows
through `build_frt_rows(meshCoor, nsmp, ...)` (`library_output.py:119`), so
columns 0-2 are physical coordinates FETCHED BY a node number and no node
number is ever written; `frt_canonical.py:131` dedupes on
`np.round(arr[:,:3], 6)` and `:150` lexsorts on those coordinates;
`compare.load_run`'s own docstring says "1 file (serial Python), 2 or 4 (MPI
Fortran), all the same". Fortran has had rank-local numbering all along, so
the committed references were produced under it.

**What is lost is CONSTRUCTION, not comparison**, and that distinction is the
whole value of the design pass. Today every index array is the gated serial
array under an injective relabelling, so a partition bug cannot hide. Three
failure classes replace that guarantee, and one — right coordinate, wrong
physics AWAY from the fault — is not reliably caught by a short coarse frt
gate. The answer is new direct-observation gates (partition arithmetic,
bitwise line slices, rank-local mesh bitwise-equal to serial under the
analytic map, per-run ownership allreduces), not a changed oracle. A design
that had proposed relaxing the comparison instead would have been refused.

I also checked the fourth claim the staging rests on, because it is the one
that would waste days: `meshgen.f90:538-545` allocates `globalOneDimCoorArr`
and `:586-596` SLICES `localOneDimCoorArr` out of it. The port must slice and
never re-accumulate — restarting the geometric stretch from a rank-local
origin perturbs coordinates by 1e-8..1e-7 m against `align`'s 1e-9 tolerance,
and the e2e failure would name node alignment while saying nothing about the
cause.

**Held for the owner, not decided here:** the memory return is UNMEASURED and
likely about half of the 1.65 GB/rank (Stage 3 is an explicit kill point), the
time return is at most 7.92 s against hundreds of ~130 ms steps, and a cheaper
alternative — 3D structured `(ix,iy,iz)` box selection on the current
build-serial-then-restrict design — buys the entire halo-surface argument
while KEEPING global numbering. Prediction registered before any build
(rule 4a): 8/16-rank ratios should fall from 1.56x/1.58x toward ~1.07x, and
**landing near 1.4x is a negative result to be reported as one** — this is the
fourth candidate explanation for the same residual after transport, placement
and element balance each measured and failed.

## HH. Landings 11-14, and the biggest one is TPV30

| commit | what | my own gate |
|---|---|---|
| `f274127`+`d5c85ef` | board refactor — `pathway_forward.md` 440,433 -> 208,343 bytes, narrative to `docs/BOARD_HISTORY.md`, items 38+89 merged, four sections | mechanical, not a read-through: 24/24 dates, 167/167 commit SHAs, 161/161 cells over 600 chars present VERBATIM in the archive, every row key still present |
| `e1888e7` | **TPV30: the port never gave a PML node its own elements' gravity** — 4.035089e+08 Pa -> 5.18e-14 | four cells run fresh by me; the Fortran side read line by line |
| `195fbd9` | perf ledger: a required `parallelism` discriminator, fail-closed | my own mutation; and I fixed the guard's own failure message, see below |
| `aca6979` | refactor round 2 — six port modules, +86/-86 | my own with/without run: frt BYTE-IDENTICAL on both backends, 752618 bytes / 1891 lines each |

**TPV30 is the session's real result.** `assembleGlobalKU.f90:44` seeds the
gravity body force into slots 10-12 of each node's TWELVE-slot PML block and
`:51-58` adds all twelve to a 12-dof node; the port added that term to
`groups[2]` alone, the 3-dof-node group sum, so **a PML node never received
the gravity of its own PML elements**. I read both sides myself before
accepting it. The onset is step 3 (t=0.125 s) — relative sigma-squared
1.08e-12, 8.18e-13, then 1.03e-05, an onset rather than growth — and the
excess is `dszz = 3*dsxx` with `dsxx = dsyy`, exactly `(lam+2mu)/lam` for
this material: a spurious VERTICAL strain rate, and gravity acts in z. Fault
tractions stayed at roundoff for 41 more steps, which is why the previous
hunt's t=1s/t=6s bracket never saw it. `test.tpv30` stays UNREGISTERED and
no reference was regenerated.

**A number I could not reproduce, and it turned out to be a label swap.** The
mission reported drv.a6 flips numpy 415 -> 332 and jax 423 -> 377. My own
re-derivation of the numpy cell read **377**, and my own jax cell read
**332** — the pair is correct, the labels were exchanged. Both inside the
unchanged 450 budget, so the gate was never in question, but I landed the fix
with the discrepancy recorded as unresolved rather than smoothed, and only
the second run settled it. Two independent cells, ten minutes, and the
alternative was a board row asserting a per-backend improvement in the wrong
column.

**I also fixed a guard I was gating.** My mutation of the ledger's e2e
emitter made `test_perf_parallelism_discriminator.py` exit 1 with a bare
`KeyError` traceback — fails closed, which is correct, but names nothing. One
line (`r.get('parallelism', '<MISSING>')`) and the same mutation now prints
`FAIL -- emitter run_e2e: {'fortran': ('<MISSING>', 4), ...}`. A guard that
cannot say what it found teaches the reader to re-run instead of read, which
is the same disease as the intermittent gate item 89 closed this morning.

**The stale-base rate is now 8 of 8, and the eighth would have been the worst
one.** Round 2's branch had gone eleven landings stale; its
`src/python/eqdyna/assembleGlobalKU.py` predates the TPV30 fix, so taking the
branch wholesale — or even taking "all the python files it touched" — would
have silently REVERTED the session's most important landing. Only its six
intended files were taken, each checked untouched on master since its base.

**Held work is not free, but landing it early would have cost more.** Round 2
sat finished for four hours because the TPV30 hunt was live in the same
directory. Landing a clarity refactor underneath an in-flight bug hunt is how
a fix gets attributed to the wrong change.

## II. Fresh 24h autopilot, resumed 2026-09-23 11:50 (predecessor stopped by user interrupt)

Grant, stated back: minor and patch tags on `master` of THIS repo, merge and
tag authority on `master` directly (the owner's 11:48 brief names v5.17.0 as
inside it). No major, no force-update or rewrite of an existing tag, no
publish, nothing outward-facing beyond this repo (item 86 and the v5.16.0
image backfill stay owner-held).

Orient, 11:49: origin/master `0077cfc`, primary clean, load 7.8/64.

**Worktrees 19 -> 4.** Fifteen reaped after checking each: every committed
branch was `git cherry` patch-equivalent to master (`-`), except
`worktree-agent-adc2fa…` (`11bb0f0`), which differs from landed `195fbd9` only
in the 7 lines of my own guard-message fix in
`test_perf_parallelism_discriminator.py`; `agent-a093ba…`'s five dirty files
were byte-identical to master (refactor round 3, landed as `090d2cf`).
**Force-release, recorded:** `agent-a2f85fa81e252b6b3` was locked by "claude
agent pid 2872539" -- that pid is dead; tree clean at `aca6979`, 0 unlanded
commits; removed with `remove -f -f`. Kept: `item32-dx250-refinement` and
`mira/jaxmpi-merged-2026-09-22` (owner-held, locked), `wei-price-tpv30` (live).

**The TPV30 pricing orphan is usable, not restarted.** Pid 1349478 is not an
orphan in the timing sense: its parent chain is intact --
`price_tpv30.py` (pid 1345968) -> `/usr/bin/time -v` (1349277) -> `run_e2e.py`
-> `python -m eqdyna` -- so the wall clock and peak RSS land in
`price_tpv30_python-numpy.log` when it exits, and the driver then runs the jax
cell sequentially. Pinned to cpu 38, box at 6-8/64 busy. Fortran cell already
in: **143.8 s wall, 0.37 GB peak RSS, max|diff| 0, SUCCESS** (tenancy 6/64).
The pricing worktree carries a scaffold registration of `test.tpv30` in
`testNameList.py` and `matrix.py` marked "never committed"; it will be
discarded, not landed -- TPV30 gating is the owner's.

**Correction to section EE.** Its mechanism paragraph ("the halo cost per rank
does NOT fall ... exactly the 1D-slab prediction ... the case for 3D
decomposition") is SUPERSEDED. Ledger lines 361-364 (sha `50c1277`, 16 ranks,
tpv104): spread placement, 2 ranks on each of the 8 NUMA nodes, jax-MPI
52.49 vs Fortran 63.27 ms/step = 0.83x; packed placement (6/7/3 on nodes
0/1/2), 96.02 vs 56.03 = 1.71x. Placement alone moves jax by -45% and Fortran
by +13% in the other direction. The "exchange" time in EE's per-rank arrays
was ranks waiting on a bandwidth-starved neighbour, not halo surface. EE's
ms/step numbers stand; every 8/16-rank ratio in that table was taken under
`least_loaded_cpus`'s packed default and is placement-biased. Board rows
(33, 43 HELD with refuted premise, 92) dispatched to zofia-kaminska.

## JJ. Resumed 15:20 after the 13:15 HTTP 429 that killed the conductor and its specialists

Grant, stated back unchanged: merge and tag authority on `master` of THIS repo,
minor and patch only; no major, no force-update of an existing tag, no publish,
nothing outward-facing. v5.17.0 is inside it.

Orient 15:20: origin/master `ef7196c`, primary clean (0 porcelain lines), CI
green on `ef7196c` (run 35899557128: build, unit-regression, e2e-ci-smoke — 3
jobs; 037170b's split jobs are GONE from `test.yml`, superseded, harmless). No
eqdyna process of ours alive; one foreign pytest at ~600% CPU. My own
`python3 testsys/run.py unit regression` at `ef7196c`: SUCCESS/SUCCESS, 3:09.

**Recovered, not redone.**
- Item 2 was split across TWO dead worktrees, not one: `agent-af3dad…`
  (Fortran emitter, uncommitted, +280/-6 in 4 `.f90`) and `agent-ac0195…`
  (schema/record/query/overhead, 848 lines, UNTRACKED — not in the brief I was
  handed). Both committed as WIP to branches, pushed:
  `origin/wip/profile-fortran-2026-09-23` (`2b78283`),
  `origin/wip/profile-guard-2026-09-23` (`e7e572a`). **Not committed:**
  af3dad's 1 ledger row + snapshot — an append (1+/0-), but stamped `6143beb`
  while produced by a MODIFIED binary; landing it would misattribute a
  profile-build timing to a SHA that never contained the code.
- Item 1 (`1102be7`, rules 24/15b/16/7) is STALE against master in at least
  three places: rule 16 lists 9 jobs / 25 of 30 cells (test.yml now has 3), and
  two "NOT YET LANDED" markers name work that landed as `6143beb` and
  `ef7196c`. Pushed as `origin/wip/rules24-2026-09-23`; back to zofia to
  rewrite against master, plus item 69's close.
- `wei-g5`'s 6 ledger rows + snapshot (tpv29/36/37 x numpy/jax at the 5 s gate
  term on clean `ef7196c`, 6/6 ok) are genuine evidence: landed `95d4220`
  (6+/0-, docs-only, rule 20b).

**Dispatched (2 concurrent, down from 3 at the 429):** zofia (item 1 + item
69), mira (item 2 emitter, all backends; pinned NUMA 2-3). The guard/A-B half
of item 2 goes to iris AFTER mira returns — the guard must be tested against
real emitter output, and a third agent in flight is what the 429 cost.
Refactor round 2's retry touches `driver.py`/`eqdyna3d.py`, which mira's
emitter will also touch: it waits for mira's landing (collision named, not
guessed).

**Mine, detached, NUMA 4-5:** TPV30 at 5 s (3 cells, 5 s reference generated
from Fortran in a scratch worktree, NEVER committed), the 5 s mutation test
(e1888e7 reverted, numpy + jax), and a pinned re-price at 20 s. The 11:37 20 s
price (Fortran 143.8 s / numpy 1020.0 s / jax 221.1 s, jax UNPINNED at 553%
CPU) is kept as the first datum.

## KK. Landings 15:30-16:05, and TPV30 priced at both terms with the 5 s mutation answered

| commit | what | my own gate |
|---|---|---|
| `d8fd06d`,`5f0131c` | rules 24/15b/16/7 (zofia): the stale 12:24 text reconciled with what `ef7196c`/`6143beb` actually landed | no `NOT YET LANDED` left for landed work; `unit regression` green on the merged tree (pinned NUMA 6-7) |
| `c33e513` | item 69 CLOSED (zofia) | row's command `test_ci_workflow_coverage.py` can go both ways (rule 14a) |
| `e93e8f7` | `test_rsfNucleation_tpv2802_td.py` (iris) — refactor round 2's retry precondition | green; aca6979's `faulting.py` -> RED `TypeError ... 9 positional arguments but 10 were given` on numpy+jax; MY OWN value mutation (drop `VINI_Z` from `backSliprate`, `faulting.py:362`) -> RED `got 0.10837790238705954, want 0.10806261170722385`; restored, sha256 `7329b52b…0f5eb1` both sides |

**TPV30, three cells per term, run sequentially by me, pinned NUMA 4-5 (cpus
32-47, verified in `/proc/<pid>/status`), box load 13-33 (foreign magic.exe x3 +
a pytest): contended, not quiet.** Scaffold registration in a scratch worktree
only (`testNameList.py`, `matrix.py`: tpv29's bound 1e-10, `abs-max`); NOTHING
committed; `test.tpv30` stays unregistered.

| term | fortran (4 ranks) | python-numpy | python-jax | 3-cell total |
|---|---|---|---|---|
| 5 s (gate) | 41.2 s, 0.39 GB, max\|diff\| 0 (self) | 266.2 s, 2.84 GB, 1.645566e-14 | 57.3 s, 4.36 GB, 1.909216e-14 | **364.7 s (6.1 min)** |
| 20 s (full) | 143.6 s, 0.39 GB, 0 | 1034.8 s, 2.84 GB, 5.184742e-14 | 177.6 s, 4.51 GB, 6.300516e-14 | **1356.0 s (22.6 min)** |

The 20 s row reproduces the 11:37 datum (143.8 / 1020.0 / 221.1 s; jax was
unpinned at 553% then, 16 cpus now). The 5 s python cells compare against a
5 s reference generated from THIS Fortran run (`frt_canonical`, 3321 rows,
sha256 `7419602702cc…98175`), scratch only — per rule 7 a committed one is its
own reviewed change, and only if the owner gates the case.

**Mutation at 5 s, answered: the 5 s gate CATCHES the TPV30 bug.** With
`e1888e7`'s `assembleGlobalKU.py` hunk reverted: python-numpy FAIL
`max|diff|=1.202240e+08 bound=1.0e-10 at row 33 col 11 (ref -1.332131e+09 vs
run -1.452355e+09)`, python-jax FAIL, identical figure. So at the everyday term
this is NOT a later catch: the bug's onset is step 3 and by 5 s it is eight
orders of magnitude past the bound on a fault-traction column. Restored,
`git status -- src` 0 lines, then the 20 s cells ran on clean code.

Gating remains the OWNER'S: this supplies price (6.1 min everyday, 22.6 min
release, peak 4.5 GB) and sensitivity (catches the one known real defect at
5 s), nothing else. Evidence:
`docs/evidence/perf-2026-09-23/wei-tpv30_price_both_terms_and_5s_mutation.log.gz`.
The scratch tree's 9 ledger rows are NOT landed: they are stamped `95d4220`
while two came from mutated source and all from a scaffold-registered tree.
That is the THIRD misattributed-SHA emission today (af3dad's, mira's 7 rows
stamped `b75c69a` from an uncommitted emitter, and these) — the ledger has no
tree-dirty field. Fix routed to iris inside item 2's record work.

**Stale comment, noted for the owner-held TPV30 decision:** `testNameList.py:8-17`
still says TPV30's divergence is "NOT yet root-caused"; `e1888e7` root-caused
and fixed it.

## LL. My own check of the profile emitter, and the first thing its placement field caught

Branch `mira/profile-emitter-2026-09-23` rebased on master (`792d0a6` in my
scratch tree), `run_e2e.py --cases test.tpv8 --backends
fortran,python-numpy,python-jax,python-jax-mpi` with `EQDYNA_PROFILE=1` then
`=0`, launched under `numactl --cpunodebind=4,5`, load ~31:
- **frt byte-identical ON vs OFF on all four backends**, 8 files (fortran
  961+961 lines; numpy 1891; jax 1891; jax-mpi 132+829+806+124), sha256 equal
  pairwise. All 4 cells SUCCESS both ways, max|diff| unchanged
  (3.051760e-11 / 3.861189e-10 / 1.983643e-10 / 1.220703e-10).
- ON: 10 profile files (4+1+1+4), all `profile_schema.validate` VALID, worst
  unaccounted/total 3.3% (python-jax serial). OFF: 0 files.
- **Defect sent back: python `fault` is 0.0 on every python backend**, folded
  into `element` (`docs/run_profile.md:69-78`). For JAX that is a real
  conflict between three owner requirements — the same buckets per rank, no new
  sync, parity absolute (a sampled un-fused step could round differently) —
  and is held for the owner as that conflict, not decided here. For NUMPY there
  is no conflict: numpy is synchronous, a `perf_counter` around the faulting
  call costs no sync, so the fold is not justified there. Re-dispatch waits for
  kai's refactor retry to leave `driver.py` (collision named).
- **Placement: `numactl` on `run_e2e.py` does NOT pin its MPI cells.** Every
  4-rank cell, Fortran AND jax-mpi, reported `cpus_allowed` rank0 [0-7], rank1
  [8-15], rank2 [16-23], rank3 [24-31] — Open MPI mapped by NUMA node over the
  parent's mask — while the serial python cells did land on 32-47. So every
  "pinned to NUMA X" statement I or a specialist made today about a 4-rank
  cell is false, including TPV30's Fortran cells in section KK (their
  wall-clock stands, their placement was nodes 0-3 including the slow cpus
  0-1) and mira's own "all pinned to 2-3". This is exactly the defect the
  placement field was built to expose, found on its first run. It could not
  be relayed to iris mid-flight (no messaging to running agents this session);
  her 4-rank A/B arms will be re-checked against their own `cpus_allowed`.

## MM. Item 2's zero-cost gate tested nothing per step; held, not landed (16:55)

Branches, none on master: `mira/profile-emitter-2026-09-23` @ `d2f83fe` (numpy
`fault` now split: tpv8 numpy `fault=0.122 s` of `element 61.2 s`, frt sha256
`d00b743f2ad6…` unchanged); `iris/profile-guard-2026-09-23` @ `c527ce2` (guard
9 checks, red both ways on real fixtures from all four backends; append-only
`docs/run_profiles.jsonl` with `tree_dirty`; e2e + run_perf wired; run_scaling
/ run_mpi_scaling NOT wired); `kai/refactor2-retry-2026-09-23` @ `bdc25fe`
(+84/-84, Td kept; drv.a6 python-jax flips 332/450 SUCCESS — same 332 as
master's jax cell this morning).

**iris's A/B (tpv8, cpus 16-31, load 32-35):** numpy ON 549.77 / OFF 518.27
ms/step, floor 17.99 -> "FAIL +31.51"; jax 100.84 / 43.18, floor 63.83 ->
"SUCCESS floor-masked"; fortran 75.47 / 73.97, floor 1.67 -> "SUCCESS +1.50";
jax-mpi 110.51 / 101.44, floor 1.75 -> "FAIL +9.07". **None of these four
verdicts measures the profiler.** `EQDYNA_PROFILE=0` gates only the end-of-run
write (`profile_emit.py:108-150`, `library_output.f90:428`); every per-step
timer the emitter added runs in BOTH arms, and per-step-by-difference cancels
the one-time write. Identical per-step code, so ON-minus-OFF is box noise —
and the two FAILs say this box's noise at load 32-35 is 6-9%, far above the
owner's 1%. iris read her FAILs as "real, needs src attention"; they are not,
and I confirmed by `grep EQDYNA_PROFILE` over `src/` before saying so.

**Decision:** item 2 does NOT land until the off switch is real (mira, in
flight: skip every ADDED timer when 0, read once at startup, OFF == the
pre-emitter step), and the A/B gate is mutation-tested the one way it has
never been: inject a known per-step cost and show the gate goes RED. A perf
gate that has only ever been observed to pass or to fail on noise has not
been shown to detect anything.

**Named conflict, for the owner:** "proven on a quiet box" vs this box today —
three foreign `magic.exe` + a pytest hold load at 30-35 all afternoon; the
OFF-vs-OFF floor alone (jax 63.83 ms on a 43 ms step) exceeds 1%. The <1%
proof may be unsatisfiable until the box quiets; the fallback I will report,
not substitute, is an operation-count bound (the emitter adds 2-6
`perf_counter`/`MPI_WTIME` calls per step, microseconds against 40-550 ms
steps) beside whatever A/B the box allows.

**Provenance flag:** `testsys/perf/profile_overhead.py:5` on iris's branch says
"DO NOT RUN THIS FROM AN AGENT SESSION (owner instruction, 2026-09-23)". No
such instruction is in the owner's 15:20 brief, which requires the A/B to
run and be gated in the perf tier. Text claiming owner authority that I cannot
trace does not land; routed back to iris to replace with the real constraint
(pinned, never concurrent with another measurement).

## NN. Speed campaign (owner top priority, relayed 17:00): the BEFORE, and the premise it breaks

**CI before:** run 35923271678 on `d2f83fe`, 4.1 min (build 0.6, smoke 1.8,
unit-regression 3.4 = critical path) — the coordinator's measurement.

**Everyday sweep before, measured:** `python3 testsys/run.py e2e` at `dc92983`
(clean worktree), 17:01:46-17:29:46, **WALL 1679.91 s (28.0 min)**, 31/31
SUCCESS, box load 35 -> 58 -> 39 (three foreign `magic.exe`, a pytest, and my
own specialists' short runs). Snapshot
`docs/perf_snapshots/e2e_cells_2026-09-23_172946_2202838.json`, 31 ledger rows
(31+/0-), log `docs/evidence/perf-2026-09-23/wei-everyday_sweep_before_dc92983.log.gz`.
Longest cells: drv.a6 numpy **1679.2**, tpv37 numpy 1572.8, tpv36 numpy 1495.1,
tpv1053d numpy 959.3, meng2023cb numpy 954.2, meng2023a numpy 934.9, tpv10
numpy 783.8, tpv29 numpy 648.7 s.

**The premise "wall ≈ slowest isolated cell; dropping tpv36/37 numpy takes
~18.5 min to ~5-6 min" does not hold on this box today, and the numbers say
why.** (1) The critical path is drv.a6 x python-numpy, not tpv36/37 numpy: it
ran 1679.2 s, longer than either. (2) Every cell ran 2-3x its solo time: tpv8
numpy 209.4 s in the sweep vs 66 s solo (mira, 16:05), tpv29 numpy 648.7 vs
~312. Total cell time ≈ 14,870 cell-s against ~25-30 free cores, and
`run_e2e.py`'s `cell_cost` bills every python cell as 1 core while a jax cell
measured 249% CPU (iris-X) — the budget of cores-4 = 60 ignores the 35 cores
already busy. The sweep is contention-bound, not critical-path-bound. So the
release-only move (landed `d488dae`) removes ~3,068 of ~14,870 cell-s (~20%) of
demand and cannot, by itself, get below drv.a6 numpy's own time; longest-first
(`e6b301d`) makes drv.a6 numpy start first, which is necessary and not
sufficient. Prediction registered before the AFTER run (rule 4a): wall ≈ drv.a6
numpy's in-sweep time, 1300-1650 s, i.e. <= 25% gain.

**Levers re-ordered by this measurement:** numpy per-cell cost across the
board (drv.a6 first — it is the floor now — then tpv36/37, meng, tpv1053d);
then tenancy-aware concurrency (bill jax at its measured cores, cap at FREE
cores) so cells stop inflating each other. mira's wedge-path profile (in
flight) is the first; drv.a6 numpy is added to the queue behind it.

**Landed:** `d488dae` (tpv36/37 x numpy release-only: SUPPORTED, not
UNSUPPORTED; everyday 29 cells, release 31) and `e6b301d` (longest-first from
the ledger's latest wall per cell; unmeasured cells first, printed). My own
gate: guard green; MY mutation (a third cell, tpv29 numpy, added to
`RELEASE_ONLY` in the file) -> `FAIL test_sweep_speed_2026_09_23 (3 check(s))`,
restored sha256 `1108e7c426e3…` both sides; `unit regression` green.

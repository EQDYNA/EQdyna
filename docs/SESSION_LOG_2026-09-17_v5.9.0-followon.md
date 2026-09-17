# Session log — 2026-09-17, v5.9.0 follow-on (wei-lin, autonomous)

Board: `pathway_forward.md`. Rules: `PROJECT_RULES.md` (17 rules). Starting
point: v5.9.0 confirmed green, tree clean. Standing order: rolling 24h budget,
report at boundaries.

## Landings this session (commit / mission / gate)

| commit | what | gate |
|---|---|---|
| `ca83335` | Docker.guide.md de-pinned `:v5.8.2` -> `:latest` (4 spots); new `testsys/regression/test_docker_guide_no_pinned_version.py` guard | fresh unit+regression green; mutation-tested (reintroducing a pin fails it); CI green |
| `a7f77d3` | item 39: `meng2023a`/`meng2023cb`/`tpv29` `CASE_BOUND` tightened `THRESHOLD` -> `1e-6`/`1e-6`/`1e-10` (stale justification "no python run ever measured" -- dead since haruto's CI_CELLS pass measured all three); item 40 opened (tpv29 numpy-vs-jax ~4.7-6x wall-clock gap, framed port-perf not CI-scheduling); item 19(c) updated (rule 17 step 3: no `TPV==34/35` branch exists in `src/fortran/*.f90`, confirmed by fresh grep) | fresh e2e re-run by wei-lin (not reused from the pre-existing comment): meng2023a/cb 4/4 SUCCESS, tpv29 2/2 SUCCESS, all values match to the last digit; unit+regression green; CI green |
| `3e13715` | tpv37 gated on all 3 backends (`mira-volkov`) — pure registration, both `faulting.f90:405`/`faulting.py:109` already had `TPV==37` branches from tpv36's own landing | mira's worktree branched from `ca83335`, one commit BEHIND `a7f77d3` (item 39's bounds) — landed via `git cherry-pick --no-commit`, not a file copy, specifically to avoid reverting item 39; diffed both shared files (`testsys/matrix.py`, `testNameList.py`) against current master after and confirmed pure additions. Fresh 6-cell oracle re-run by wei-lin (tpv37+tpv36, 3 backends): all 6 match her report to the last digit. Unit+regression green; CI green. Worktree reaped clean (`e4d9f91`, no uncommitted work) |

## Model-override provenance (item 24, tpv30 unblock)

`dunyu-liu` dispatched on the default Fable-5 model **four consecutive times**
for item 24(b)(c)(e)(f), each an immediate `429` ("You've reached your Fable 5
limit"), spaced ~30-90 min apart. Tested and disproved a concurrency
hypothesis: retried alone after the only other in-flight agent (mira) finished
— still 429'd immediately. This pointed at a Fable-5 session/account budget,
not a per-request or concurrency limiter.

**Per the coordinator's instruction, re-dispatched `dunyu-liu` with an explicit
`model: "opus"` override** (opus preferred over sonnet per the coordinator —
item 24 is genuine numerical/config work on the viscoplastic path, not a
mechanical change). This went through cleanly (worktree created, agent
running, no 429). **Any commit landing from this mission was authored under
an opus override of dunyu-liu's normal Fable-5 default — recording this here
so it does not vanish from the provenance trail; the landing commit message
must also say so.**

## Deliberately not dispatched this session

- **tpv34/tpv35**: rule 17 step 3 done (no `TPV==34/35` Fortran branch exists
  — net-new physics on both sides, not a port). Step 1 (spec fetch) is cheap
  and can run on whatever model budget is available next; Python scoping must
  not start before the Fortran branch question is settled.
- **Item 33** (NUMA re-measurement): box was measured genuinely idle (load
  1.4-1.9) mid-session, but my own concurrent oracle re-runs and agent
  dispatches were themselves contending for it the whole time — never a
  clean, uncontended window this session. Do not measure NUMA scaling against
  your own load.

## Owner decisions in force, not reopened

Item 29 closed (stashes dropped, confirmed empty `git stash list` this
session). Item 36 decided: no history rewrite. Item 17 deferred, items 7/9/10
behind it. Working TPV joins the default sweep. All backends or it does not
land. tpv29 50m closed.

# Session log — 2026-09-18 autopilot (v5.11.1)

## Budget read

Board nearly exhausted per prior session's own closing note ("no other board
items are actionable outside perf... and the tpv30/drv.a6 hands-off zone").
One substantive target: item 33's remaining question (why jax falls off past
16 cores). Three hygiene rows (16, 35, 38) named as fallback-only, not to be
chased unless item 33 stalls. Grant: patch/minor on non-default branch per
default autonomous-mode reading; nothing in this run's scope needs wider than
that (no release cut planned unless item 33 lands a real fix).

## Orientation

HEAD `cacf738` (master), clean tree, matches v5.11.1. Box: load ~quiet
(mpstat: NUMA nodes 4-7 ~99-100% idle; only 2 of the owner's `train.py`
processes actually running at ~1 core each, not 4 as briefed — doesn't change
the "box is free" conclusion). Confirmed live trap: ambient shell
`PATH`/`VIRTUAL_ENV` resolves to the owner's unrelated `gns` venv
(`/home/utig5/dliu/gns/gns/venv_cotopaxi`) whose jax defaults to CUDA GPU
devices. No dedicated EQdyna venv exists; the project's own perf tooling
(`testsys/perf/run_scaling.py:322`, `run_numa_scaling.py:222`) already forces
`JAX_PLATFORMS=cpu` on that same venv rather than using a separate one — used
that exact invocation for the dispatch below rather than guessing a new one.

Checked item 35 (duplicate same-second CI runs) for recurrence out of
curiosity while the dispatch was in flight: `gh run list --limit 30` shows 4
more shas with same-second duplicate push-triggered runs since the item was
last touched (`f607bfd8`, `e588708e`, `d47876b4`, `a5c6244f`). Confirms the
item's own prediction ("cheap to reproduce next time it happens"). Not
investigated further per that row's existing disposition — cost is wasted
runner-minutes only, no incorrectness. Left for whoever picks up hygiene
rows; not chased this session.

## Item 33 dispatch — jax CPU scaling plateau, op-level diagnosis

Dispatched dunyu-liu (worktree-isolated) to build a controlled microbenchmark
distinguishing scatter-add-specific vs memory-bandwidth-general causes of the
2.52x-at-16/1.98x-at-32 plateau, per the mission brief in the prompt.

First dispatch (agentId a7288ece, default model) failed before creating a
worktree: `429, model claude-fable-5, rate_limit`. Per this repo's own
recorded trap ("A 429 is model-specific — override the model, record which
ran what"), re-dispatched immediately with `model: sonnet` (agentId
a437a272). No partial work from the failed attempt to salvage (it died before
`git worktree add`, confirmed via `git worktree list`).

Agent completed (77949 tokens, 22 tool calls, 454.7s). Report: verdict "(a),
refined by HLO evidence" -- scatter-add specifically fails to CPU-parallelize
(0.75x-0.85x speedup 8->32c across 2-3 runs) while a matched-bytes control op
"never regresses" (1.03x-1.76x); HLO partition counts identical for both ops
(3@8c, 6@32c) so framed as inter-socket write contention, not partition
starvation.

**Re-verified myself before landing (gate axis 3) -- caught an overclaim.**
Re-ran the identical committed script 3 more independent times, same sha,
same interpreter: one run gave control speedup = 0.74x, directly falsifying
"control never regresses." Applying the script's own verdict threshold
(control>=1.5x AND scatter<0.6*control for verdict a) to all 5 completed runs
(2 original + 3 mine): only the FIRST run meets that bar; the other 4 read as
"(c) inconclusive" by the identical rule, and one of mine (control=0.99x,
scatter=1.02x) inverts the claimed asymmetry. The box's own busy-check
printed a live multi-tenant roster each run (cisco-amp-scan-svc, cosmo,
cshyu, enze, junlin, lp, markw, messagebus, rdomeyko, ron, yukoo, zjia) -- a
genuinely shared login node, not the "finally quiet" box the dispatch brief
assumed (which named only the owner's 2 train.py jobs).

Corrected `NOTES_item33_scatter_probe.md` in place (marked and dated, original
text preserved for the record, not silently overwritten) to: **NOT SETTLED**.
The HLO partition-count finding is a compile-time property, unaffected by
runtime noise, and stands as-is. The timing asymmetry is at most a weak,
directionally-suggestive signal on a box whose noise floor is comparable to
the effect claimed.

Landed `036b4ba` on master (`testsys/perf/probe_scatter_bandwidth.py`,
`NOTES_item33_scatter_probe.md`, `testsys/perf/scatter_bandwidth_last.json`)
-- new standalone files, base = current HEAD confirmed via `git diff --stat
cacf738`, zero staleness. Direct master landing (small, report-only, no
gate/kernel surface); noted here as the required non-default-branch departure
disclosure per this session's contract, matching this project's own
established convention of direct-to-master patch landings all session.

Fast gate on the landed commit: `python3 testsys/run.py unit regression`,
SUCCESS both tiers, exit 0.

Dispatched zofia-kaminska (rule 19 -- I do not author board rows) with my
exact verified numbers to append the item 33 row update. Verified her diff
myself before landing: `git diff --stat pathway_forward.md` = 1 file, 1
line changed -- only the item-33 row's own line, no other row or file
touched, content matches what I commissioned. Landed `90fdc12`.

## Item 33 status at end of session

NOT CLOSED, sharpened. Standing gap unchanged in kind but now more precise:
an idle/reserved box is still needed, AND the next probe should average
multiple repeated trials per configuration (this one did single-shot timing
per config) alongside the still-unrun `perf c2c`/coherence-counter check
named in the probe's own open-questions section, before either confirming or
killing the write-contention hypothesis.

## Close

Board confirmed exhausted for this run: item 33 is the only substantive
target and it is now sharpened rather than closed; the remaining P2-P4 board
rows are CLOSED, DEFERRED-BY-OWNER-DECISION, or explicitly hands-off
(tpv30/drv.a6). Hygiene items 16/35/38 remain fallback-only per the mission
brief and were not chased beyond the incidental item-35-recurrence check
logged above. No fabricated work. Stopping on genuine queue exhaustion, per
the standing instruction to do exactly that rather than manufacture a task.

Tags this session: none (no version bump -- this was tooling/diagnosis, not
a release point; matches this project's own precedent of batching several
small commits between patch tags, e.g. `9b0e212`/`c6c31db`/`cacf738` all
landed after v5.11.1 with no intermediate tag).

Worktrees: `agent-a437a272308248b00` (item 33 mission) -- fully committed and
cherry-picked to master (`2ad69c4` -> `036b4ba`), safe to reap, no unlanded
work in it. Reaping now.

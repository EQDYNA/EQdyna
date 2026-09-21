# Session log — 2026-09-21 — jax 16-core scaling (pathway item 43)

Conductor: wei-lin. Start: 11:01. Base: `51b0fa1` = `v5.13.0`, tree clean,
`origin/master` level (verified, not assumed: `git status --porcelain | wc -l`
= 0; `git log --oneline -1 origin/master` = `51b0fa1`).

ONE objective for 24h: jax backend scales ~linearly to 16 cores (bar 12-14x),
against Fortran's measured 14.39x. Closes item 43. Not a board sweep.

## Budget reading (stated before spending any)

- **Wall clock is the scarce resource, and it is mostly Mira's, not mine.** The
  work is one deep mission (shard_map over a 16-device CPU mesh) already in
  flight. My spend goes to the merge boundary and to keeping the loop alive,
  not to re-deriving her measurements.
- **Cores are the second scarce resource and they are contended.** `nproc` 64,
  but load average at start was 19.45 / 42.03 / 45.54 from a foreign tenant
  (user `kehua`, ~21 single-core jobs — the same contamination items 33/42/43
  all record). Every scaling number this campaign produces must be read against
  a recorded contention sample; a number without one is not evidence.
- **Therefore heavy jobs SERIALIZE.** Item 42: three timeout near-misses across
  two sessions, all from self-inflicted stacking. Mira's measurement passes and
  my gate sweep (`run.py all`, 30 cells, 36-90 min) may never overlap. I hold
  my sweep until she reports a landing candidate.
- **Tag authority as granted, restated for correction before the first tag:**
  patch and minor on `master` of this repo. No major, no publish, no
  force-update of an existing tag. Campaign sign-off is the owner's.

## Agents

- `mira-volkov` `ad64043fb3e42e07a` — ADOPTED, not respawned. Worktree
  `/home/utig5/dliu/wt-jaxshard`, branch `jax-shardmap`, created 10:59
  (2 min before my first probe), base `8919405` = `v5.13.0`.
- **Base staleness assessed, not waved past (gate axis 4):** her branch is 2
  commits behind `51b0fa1`, and both are DOCS ONLY — `ec4d464`
  (`docs/SESSION_LOG_2026-09-19_autopilot.md`, +88) and `51b0fa1`
  (`pathway_forward.md`, +4). Zero source overlap, so her code work is
  unaffected. Consequence recorded here because it is the one way this bites:
  a wholesale copy of `pathway_forward.md` out of her worktree would revert
  those 4 lines. Board rows are zofia-kaminska's (rule 19) and no mission
  writes that file, so the exposure is closed by ownership, not by luck.
- **Channel limitation, stated plainly:** this deployment has no
  `SendMessage`/`ListAgents` tool — only `Read`/`Edit`/`Write`/`Bash`/`Agent`.
  I cannot query her for status and cannot enumerate live agents. Adoption is
  therefore by observation of artifacts she owns (worktree mtimes, git index,
  commits, `NOTES_*`, her processes), which is the liveness definition I use
  anyway. This is the one thing about this run I would ask the owner to fix:
  without a message channel I cannot relay a mid-flight correction to her.

## Gate shape for this objective (decided up front)

A `shard_map` rewrite lands in `src/python/eqdyna/backend.py`, whose 3 scatter
sites are reached from 8 call sites in `assembleGlobalKU.py` and 1 in
`faulting.py`. That is every python cell in the table, both backends. So:

- Merge gate = full `python3 testsys/run.py all` (30 cells), unbounded or
  >=90 min, nothing else heavy on the box. Not a subset.
- Parity first: every cell inside its own `testsys/matrix.py` bound. No
  reference regenerated (rule 7), no bound loosened.
- jax-CPU only. No gate may assume GPU bit-identity (8.9e-08 run-to-run).
- Per-step cost by DIFFERENCE over two step counts, both backends — the bias
  `88f5227` removed (it read 91.25 vs a true 64.74 ms/step at 16 ranks).
- Every verdict carries a property of the compared content (cell count, row
  count, a grep count). "Green" with no such property is a vacuous gate; three
  of those passed in the last campaign.

## Stopping rule (owner's, recorded so it is not softened later)

Stop on a named mechanism WITH a measurement: "plateaus at Nx because <thing>
costs <measured amount> per step." The withdrawn 6x kill criterion does not
apply. The Amdahl 6.7-10x cap in item 43 was computed off AUTO-partition counts
and does not bind explicit sharding — it is not grounds to decline the attempt.

## Timeline

- **11:01** — Orient. Git state verified clean and level. Item 43 read in full.
  Mira confirmed live (worktree 2 min old, git index advancing) — not
  respawned. Contention baseline recorded above.

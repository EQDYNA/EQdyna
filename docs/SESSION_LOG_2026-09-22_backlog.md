# Session log — 2026-09-22, backlog clearing (wei-lin)

Mandate: clear the open board items that do not need the owner. Owner-gated
rows (65, 54, 63, 57, 32, 17, 19(b)) were not touched; nothing was tagged,
nothing was published, `VERSION` still reads 5.15.0 and is still untagged.

Board rows are `zofia-kaminska`'s (rule 19). This log records what landed and
what gated it; the rows were handed to her with the literal command output.

## Landings, in order

| time | commit | what | gate |
|------|--------|------|------|
| 13:5x | `25453e5` | item 50 — `testsys/regression/test_pretag_ci_negative.py`, all six pre-tag outcomes; CI `unit-regression` gains `fetch-depth: 0` | `run.py unit regression` exit 0; new test exit 0, 8 cases; mutation-checked RED |
| 14:0x | `667dba0` | item 49 — `docs/HISTORY.md` + `test_history_table.py` (generator and verifier) | `run.py unit regression` exit 0 |
| 14:1x | `1b896fc` | fix for `667dba0` — the generated doc re-derived the repo's own size | `run.py unit regression` exit 0 |
| 17:2x | `d4f625e` | item 61 — per-rank XLA compilation-cache directory (cherry-pick of `14f0b7f`, code+tests only) | full `run.py all`: 31/31, 9413.8 s |
| 17:2x | `2c5c88c` | evidence for `d4f625e` + 33 ledger rows | `test_perf_ledger.py` exit 0 |
| 19:1x | `a16cc98` | item 58 — one GPU cell per card, preallocation stated not inherited | full `run.py all`: 31/31, 3340.2 s; GPU before/after |
| 19:2x | `1d8c70d` | evidence for `a16cc98` + 33 ledger rows | `test_perf_ledger.py` exit 0 |

Re-verified against their own evidence commands, no code change needed
(the board already carried them as resolved and the code agreed): item 23
(`test_drucker_prager_kernel.py` exit 0), item 51 (`grep -c` nonzero in both
files), item 52 (`test_equilibrium_dump.py` exit 0).

## Deviations from the brief, recorded rather than asked about

**Item 23 was NOT routed to `iris-vermeulen`.** The brief asked for the
isolated Drucker-Prager kernel test to be designed; the board already carries
it as resolved at `d5da43e` and the test exists and passes. Dispatching would
have re-bought work already done. The plan was stale; the code was
unambiguous.

**Item 61's diff touches two files the brief fenced off** (`backend.py`,
`driver.py`), because the brief's "everything below avoids that surface" is
not true of `14f0b7f` — it is precisely a `backend.py`/`driver.py` change.
The explicit instruction to land it wins over the general fence, and landing
it strictly reduces risk to the peer measuring on those files, since the
wedge it removes is triggered by concurrent jax. Flagged upward.

**Items 50, 51, 52 were already VERIFIED/LANDED on the board** when the brief
described them as open. Re-ran their evidence rather than re-doing the work.
For item 50 that re-run is what surfaced the real finding: two of the six
outcomes had stopped reproducing, which is why the fix is a committed test
rather than another by-hand campaign.

## Contradictions and bounds worth carrying forward

**Item 50 — a "verified" row can decay.** The 2026-09-21 campaign drove all
six outcomes from live CI history, correctly. Today, `--pre-tag 154aa0f`
returns PASS, not PENDING: the run it used has completed. A verification whose
evidence command depends on a transient remote state is not a standing gate.

**Item 58 — the mechanism is not "two cells cannot share a card".** With two
concurrent processes on a 40 GB A100 the OOM does NOT reproduce; XLA's 0.75 is
of FREE memory and the second still fits. It took FOUR to reproduce it. The
fix is right either way, but the sentence "two concurrent cells cannot both
hold 75%" is not what the hardware does.

**Item 61 — the hang was not reproduced.** It is intermittent and has no
committed repro. The deterministic coverage is three unit tests, not the
sweep. Anyone reading the green sweep as proof that the wedge is gone is
reading more than it says.

**Tenancy.** `gate-d4f625e` ran with several other agents' jobs on the box
(load ~100 on 64 cores) and took 9413.8 s for the same 31 cells that took
2904.6 s on a quiet box and 3340.2 s in the second sweep. Parity verdicts are
load-independent; the wall clocks in those ledger rows are not comparable
across tenancy. Both evidence READMEs say so.

## Process notes

Both full sweeps ran DETACHED in an isolated worktree
(`.claude/worktrees/wei-item61`, removed at the end), polled by artifact, not
in the main checkout — which turned out to matter: two other agents were
committing to the main checkout during this session, and HEAD moved under a
command mid-way through. Three worktrees existed that the brief did not
mention.

One self-inflicted failure worth the entry: the first idle-detector for
"is the box busy" used `ps ... | grep <pattern>`, which matched the polling
shell's OWN command line, so a gate waited 8.5 minutes for a box that was
already free. `pgrep -x <comm>` cannot match the poller. Same family as the
`pkill -f` lesson.

# Gate evidence — v5.16.0

## The gate that passed

`run_all_sweep_2026-09-22.log.gz` — the full `python3 testsys/run.py all` that
gates v5.16.0.

| | |
|---|---|
| SHA gated | `88a4012` |
| where | an ISOLATED worktree, `.claude/worktrees/wei-gate-v5160`, detached at that SHA |
| build | `./install-eqdyna.sh -m ubuntu`, exit 0 |
| start / end | 2026-09-22T22:25:50-05:00 → 22:55:51-05:00 (1801 s wall) |
| unit | SUCCESS (exit 0) |
| regression | SUCCESS (exit 0) |
| e2e | SUCCESS (exit 0) — **31 of 40 cells in the 10 case x 4 backend table, 31 passed, 0 failed** |
| ledger | 31 cell-wall-clock rows appended; snapshot `docs/perf_snapshots/e2e_cells_2026-09-22_225550_2502752.json` |
| exit | 0 |

`88a4012..<release commit>` touches only `PROJECT_RULES.md`,
`pathway_forward.md`, `docs/` and the release files themselves — no solver
source, no reference, no `testsys/` file. The gated code surface is the tagged
code surface.

## The gate that did NOT pass, kept deliberately (rule 8)

`collided_sweep_2026-09-22_INCIDENT.log.gz` — the FIRST attempt at this same
gate, run in the MAIN CHECKOUT starting 21:48, and destroyed by a concurrent
session. It is preserved because it is the worked example behind
`PROJECT_RULES.md` rule 21a, and because the shape of its failure is the point.

What happened: `testsys/e2e/run_e2e.py:587-594` rotates `$REPO_ROOT/test` to
`test.prev` at startup, unconditionally and with no lock. At 22:12:20 a second
Claude session ran its own e2e in the same checkout (its artifact,
`docs/perf_snapshots/e2e_cells_2026-09-22_221220_2462614.json`, landed on
master in `88a4012`) and rotated this sweep's in-flight tree out from under it.

The damage, as it appeared in the log:

    FileNotFoundError: [Errno 2] No such file or directory:
      '/home/utig5/dliu/EQdyna/test/test.tpv1053d.python-numpy/frt.txt0'
    -- cell: test.tpv1053d x python-numpy (1500.1s) --
    FAIL test.tpv1053d x python-numpy

raised at `src/python/eqdyna/library_output.py:120` from `eqdyna3d.py:522`,
which writes by ABSOLUTE path — so the rename broke it at the very end of a
25-minute run. Twenty-five cells had already passed; four more
(`test.tpv29`, `test.tpv36`, `test.tpv37`, `test.drv.a6`, all python-numpy)
were still running into paths that no longer resolved and were doomed for the
same reason. They were killed by PID.

**The reason this is kept: it is a FALSE RED that is indistinguishable at a
glance from a real solver failure.** `FAIL test.tpv1053d x python-numpy` in a
sweep summary is exactly what a genuine parity regression looks like. The
proof that it was not one is in the passing log beside it — the same cell,
same SHA, run in isolation twenty minutes later:

    -- cell: test.tpv1053d x python-numpy (907.0s) --
       max|diff|=4.569608e-06 bound=1.0e-04 at row 1664 col 13, 4005 fault nodes compared
    SUCCESS test.tpv1053d x python-numpy

A conductor who had taken the first result at face value would have reverted a
good change or re-gated a landing that was already fine.

The violation was mine: rule 21a did not exist yet, but the isolation
discipline it codifies did, and I ran a 40-minute gate in the shared main
checkout anyway. The mechanical fix — a lock on `$REPO_ROOT/test` — is
pathway item 70 and is NOT done; until it lands, rule 21a is enforced by
reading it.

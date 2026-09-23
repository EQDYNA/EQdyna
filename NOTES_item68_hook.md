# item 68 — pre-commit guard for rules 21 / 21b (worktree checkpoint log)

Worktree: `/home/utig5/dliu/EQdyna/.claude/worktrees/iris-item68`
Branch:   `iris/item68-precommit-hook`, branched from 894cdc1 (v5.16.0).
Scope:    hook + installer config line + regression guard. No solver source,
          no reference, no bound touched. Not merged, not pushed.

## Step 1 — `testsys/hooks/pre-commit` (new)

Refusal condition, and nothing else (rule 21b's amendment settles this — no
env marker, no role test):

    realpath(git rev-parse --git-dir) == realpath(git rev-parse --git-common-dir)

Both are resolved with `cd ... && pwd -P` before comparing, because git returns
the bare string `.git` for both at the top of the main checkout but absolute
paths from a subdirectory or a linked worktree; comparing raw strings would be
right by accident in one case and wrong in the other. If either `rev-parse`
fails, or either path fails to resolve, the hook REFUSES rather than assuming
it is in a worktree (rule 2 — no fallback).

Verified in a throwaway sandbox: refuses in the main checkout (exit 1, HEAD
unmoved), allows a linked-worktree commit (exit 0), and `git merge --ff-only`
still works in the main checkout because a true fast-forward creates no commit
and never invokes pre-commit.

## Step 2 — `install-eqdyna.sh`

One block added immediately before the rule-13 chmod block, inside the
existing `if [ -n "$MACH" ]` body so it runs for `-m` and `-c` alike:

    git config core.hooksPath testsys/hooks

guarded by `git rev-parse --git-dir` so a non-git source tree says so loudly
instead of silently skipping. The value is RELATIVE on purpose: git runs hooks
from the top level of the working tree, so each linked worktree uses its own
checked-out copy of `testsys/hooks/` — one install covers every worktree,
present and future. Nothing else in the script was touched; exec bit unchanged
(index mode still 100755, cf. c4b13e5).

Idempotency measured, not asserted: two consecutive `bash install-eqdyna.sh -c
ubuntu` runs in a sandbox repo left `.git/config` byte-identical
(md5 0217caa5f3ddf54e365e7fb86978872e both times) with exactly one
`core.hooksPath` entry.

## Step 3 — `testsys/regression/test_precommit_main_checkout_guard.py` (new)

Asserts, in the existing regression style (standalone script, own SUCCESS/FAIL
banner, non-zero exit, no skips):

1. the hook is TRACKED (`git ls-files -s`, mode read from the index);
2. it is EXECUTABLE — index mode 100755 AND the working-tree bit (git silently
   ignores a non-executable hook, which is the failure that reads as installed);
3. it ACTUALLY REFUSES a main-checkout commit, with "REFUSED" and a worktree
   instruction in the message, and HEAD does not move;
4. it ACTUALLY ALLOWS a linked-worktree commit, and leaves `merge --ff-only`
   working — the load-bearing direction, since a hook that blocks worktree
   commits stops every agent in this project;
5. `install-eqdyna.sh` sets `core.hooksPath` (a hook nobody installs is inert).

Behaviour is exercised in a throwaway `git init` under `tempfile.mkdtemp()`
with its own `core.hooksPath`. This repository's shared config is never
written — it is shared by every live worktree and other sessions are
committing through it right now.

Mutation-verified (each mutation applied, guard re-run, mutation reverted):

| mutation | guard |
|---|---|
| hook index mode 100644 + working-tree 644 | FAIL (both facets named) |
| hook body replaced by `exit 0` | FAIL (commit allowed, HEAD moved) |
| hook body replaced by unconditional refusal | FAIL (worktree commit blocked) |
| `git config core.hooksPath` line deleted from installer | FAIL |
| none (baseline) | SUCCESS |

## Step 4 — gate

`python3 testsys/run.py unit regression` in this worktree; result recorded in
the session report.

## What activating this breaks

Nothing, until someone runs `install-eqdyna.sh` in a given clone: git installs
no hooks on clone and `core.hooksPath` is unset today, so the tracked hook is
inert in every existing checkout, including the live ones. After an install,
the only action that changes is `git commit` from the main checkout, which
rules 21 and 21b already forbid. `git merge --ff-only`, `git commit` in any
worktree, `git add`, rebases inside a worktree, and CI (which clones and never
runs the installer before committing — it does not commit at all) are all
unaffected.

One residual, recorded rather than fixed: a worktree checked out on a branch
that predates this commit has no `testsys/hooks/pre-commit`, and git treats a
missing hook file as "no hook" rather than an error. Such a worktree is
unguarded until it carries the file — which is harmless (worktree commits are
permitted anyway) and only matters for the main checkout itself, which does
carry it once this lands on master.

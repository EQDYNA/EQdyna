# Item 75 (P3) -- rule 21c made mechanical: one writer for the board, enforced as SEPARATION

Branch `item75-one-writer-hook`, from `51f20e0`. Worktree only; not merged, not pushed.
`PROJECT_RULES.md` and `pathway_forward.md` are NOT touched by this change -- editing them
here would be the exact act the guard exists to refuse. Proposed row text is in the report.

## What landed

1. `testsys/hooks/pre-commit` -- ONE hook, extended, not a second one.
   The main-checkout refusal (item 68, rules 21/21b) stays FIRST and its message is
   byte-identical; only its `if` was inverted so the file can continue past it.
   New second check: refuse when the staged set (`git diff --cached --name-only`,
   against the empty tree when there is no HEAD) contains `PROJECT_RULES.md` or
   `pathway_forward.md` TOGETHER WITH any other path. The message names the other
   files and gives the two-commit split.

2. `testsys/regression/test_precommit_board_separation_guard.py` -- sibling guard.
   Drives real `git commit`s in a `tempfile.mkdtemp()` sandbox. This repo's shared
   `core.hooksPath` is never written (hard constraint: other live worktrees share it).

3. `testsys/check_board_separation.py` -- the merge-boundary half, added because a
   measurement forced it (see the amend gap below). Range is REQUIRED; no default.

## Measured, not assumed

| case | exit | consequence |
|---|---|---|
| `git merge --ff-only` | 0 | creates no commit, hook never runs -- unaffected |
| automatic merge commit (`git merge --no-ff` of board branch into code branch) | 0, 2 parents | **git does NOT run `pre-commit` for it** -- the merge path is NOT jammed |
| `git commit --amend` on a board-only commit | 0 | allowed |
| `git commit --amend` adding code to a board-only commit | 0 | **BLIND SPOT** |

The blind spot is structural: `git diff --cached` compares the index against the
commit being REPLACED, so an amend shows the hook a clean staged set while the commit
that lands is mixed. Git passes `pre-commit` no argument and sets no marker that
distinguishes an amend, so it cannot be closed in the hook -- closing it would require
guessing (reading `/proc/$PPID/cmdline`), which rule 2 forbids. It is closed instead
over finished history by `testsys/check_board_separation.py`, which is where the
2026-09-23 incident was actually caught.

Stated limit of the range checker: merge commits are SKIPPED (an automatic merge
introduces no content of its own; its first-parent diff is the whole other branch).
An evil merge -- a merge commit hand-edited to carry content -- is not caught.

## Mutation verification (both new checks are load-bearing)

- Hook `case` arm changed to `NOSUCHFILE_A.md|NOSUCHFILE_B.md)` -> guard FAILS with
  6 findings, including "HEAD MOVED despite the mixing refusal". Restored.
- Range checker `if board and other:` -> `if False:` -> guard FAILS with
  "check_board_separation PASSED a range containing a commit that mixes". Restored.

`python3 testsys/run.py unit regression` -> exit 0 with both guards in the tier.

## Open

- Rule 21c says AUTHORSHIP; the hook can only see SEPARATION. Suggested rewording in
  the report so the rule and the mechanism describe the same thing.
- `check_board_separation.py` is not wired into CI. Wiring it (`origin/master..HEAD`
  on PR jobs) is the remaining half and needs a conductor decision on the range.

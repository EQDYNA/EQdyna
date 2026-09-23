# Item 80 — wiring `testsys/check_board_separation.py` into CI

Branch `item80-ci-wiring`, worktree only, never pushed. Base `4ee171b`.
Filed under `docs/` per the convention `4ee171b` itself established
("docs: relocate item-74 mission notes under docs/").

## What landed

ONE step, `Check board/code commit separation (rule 21c)`, first in the
existing `unit-regression` job of `.github/workflows/test.yml`. No new job,
no new runner, no change to any other job's steps and no change to the e2e
matrix. That job is the only one of the seven already carrying
`fetch-depth: 0`, which the range needs; the step is placed ahead of the
apt/pip installs because it needs only `git` and the runner's system
`python3` (the checker is stdlib-only), so a violation costs seconds.

## Range

| event | range |
|---|---|
| `push` | `github.event.before..github.sha` |
| `pull_request` | `pull_request.base.sha..pull_request.head.sha` |
| anything else | REFUSE, exit 1 |

Expressions are interpolated into `env:`, never into the shell body
(expression substitution inside `run:` is a script-injection surface).

Two push cases `github.event.before` cannot express, neither allowed to pass
silently (rule 2):

- **first push of a branch** — `before` is 40 zeros;
- **force-push** — `before` names a commit this clone no longer has
  (detected by `git cat-file -e`, not assumed).

Both take the STATED fallback `origin/master..HEAD`, which on a feature
branch is a SUPERSET of the push — every commit the branch introduces — i.e.
wrong in the safe direction, never narrower. `fetch-depth: 0` is what makes
`origin/master` resolvable (it fetches `+refs/heads/*:refs/remotes/origin/*`).

Two ways that fallback could itself check nothing, both closed:

- `origin/master` absent → REFUSE, exit 1.
- `origin/master..HEAD` EMPTY (HEAD already an ancestor of master — the shape
  a force-push to master leaves) → narrow to the head commit alone
  (`HEAD^..HEAD`) rather than report a vacuous zero; refuse if HEAD is a root
  commit. Demonstration 7 below is the mutation proving this path still
  fails on a mixed tip, not just prints a number.

## Demonstrations (literal, dry-run)

The `run:` block was extracted from the YAML by `yaml.safe_load` and executed
verbatim — the YAML text itself was run, not a reimplementation of it —
against a synthetic repo carrying one MIXED commit (`pathway_forward.md` +
`src.f90` in one commit), clean code commits and board-only commits.

| # | scenario | range used | exit |
|---|---|---|---|
| 1 | push, normal `before`, range contains the mixed commit | `CLEAN..MIXED` | **1 FAIL** |
| 2 | push, normal `before`, clean range | `CLEAN..CLEANTIP` | 0 PASS |
| 3 | push, `before` = 40 zeros (first push) | fallback `origin/master..HEAD` | **1 FAIL** |
| 4 | push, `before` = missing object (force-push) | fallback `origin/master..HEAD` | 0 PASS |
| 5 | `pull_request`, base..head clean | `base..head` | 0 PASS |
| 6 | push on master, zeros, HEAD == origin/master, clean tip | narrowed `HEAD^..HEAD` | 0 PASS |
| 7 | **mutation of 6**: same path, tip is MIXED | narrowed `HEAD^..HEAD` | **1 FAIL** |
| 8 | fallback needed, `origin/master` absent | none | **1 REFUSED** |

## Cost, measured

- 90 ms per invocation (10 runs, 5-commit range, this repo, wall 0.898 s).
- 5.36 s on the whole 781-commit history — the absolute upper bound, never
  a range this workflow constructs.
- Against a ~20-minute CI critical path (`e2e-ci-python-tpv29`), and inside
  an existing job, this is not a meaningful addition.

## Corpus count before trusting the gate

Full history: **54** mixed commits in 781. Recent: `HEAD~10..HEAD` 0,
`HEAD~30..HEAD` 2, `HEAD~60..HEAD` 3, `HEAD~120..HEAD` 7. So the gate is
green on current practice and would be red on older history — which no range
the workflow builds ever reaches, because both the primary ranges and the
fallback cover only commits NOT yet on master.

## Stated limits

- Merge commits are skipped by the checker, deliberately. An **evil merge**
  — a merge commit hand-edited to carry a board file plus code — is NOT
  caught. This is written into the workflow comment, not just here.
- A board-ONLY push never reaches this step: the workflow's `paths-ignore`
  starts no run for one. Correct and free — a MIXED commit by definition also
  touches a non-board file, so it always triggers.

## Not done / reported instead

- `pathway_forward.md` has **no item 80** at `4ee171b`, on `origin/master`,
  or in the main checkout. Not written here (rule 21c: one writer).
- No regression guard asserting this step exists in the YAML. `testsys/
  regression/test_ci_workflow_coverage.py` guards e2e cell coverage only.
  Recommended, deferred to the owner as scope.

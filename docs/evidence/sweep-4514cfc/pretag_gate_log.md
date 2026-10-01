# v5.21.0 pre-tag gate log

Rule 15f step 1 ("board row first, before the release PR merges") was missed
in real time: the release PR (#73, M=`4514cfc`) and the sweep-evidence commit
(T=`bb91f5a`) both landed before `pathway_forward.md`'s dated Tasks-done row
for v5.21.0 existed anywhere in this branch's history. This file is the
transcript of discovering and closing that gap, kept as evidence (this
directory, `docs/evidence/`, is not `paths-ignore`'d, unlike
`pathway_forward.md` itself -- see rule 15g) rather than silently retried.

## What happened, in order

1. `python3 testsys/regression/check_pretag_ci.py --pre-tag bb91f5a8bc649f0c37d883834ed882c9ff290bd7`
   returned exit 6, `RELEASE_DOCS_NOT_READY`: CI was green and the sweep was
   sufficient, but `bb91f5a`'s own tree had no dated Tasks-done row naming
   v5.21.0.
2. The row was added to `pathway_forward.md` and pushed directly to master
   (rule 25's board carve-out) as `7718c3a`. That commit touches ONLY
   `pathway_forward.md`, which is in `test.yml`'s `paths-ignore` list (rule
   15g: it is pure prose, read by nothing under `testsys/` as physics-gate
   input), so it could never trigger a `push`-event CI run of its own.
3. `check_pretag_ci.py --pre-tag 7718c3aead73ea3cd2fca21b98a2719704738576 --ack-paths-ignored-parent`
   returned exit 3, `PATHS_IGNORED`: the 2026-09-30 hardening (the v5.20.2
   incident) deliberately removed the ancestor-evidence escape hatch this
   flag used to grant. The tool's own remedy text: "push a commit that also
   touches a non-ignored path, or narrow the workflow's paths-ignore list".
4. Narrowing `paths-ignore` (`.github/workflows/test.yml`) is a `testsys/`
   `.github/` design change gated by rule 25 and explicitly flagged by rule
   15g as a judgment call needing its own review, not a unilateral mid-release
   edit -- out of scope here.
5. This file is the "commit that also touches a non-ignored path": a true,
   non-fabricated record of this exact gate interaction, placed under
   `docs/evidence/sweep-4514cfc/` (the sweep this release's tag is justified
   by) rather than `docs/notes/` (also paths-ignored) so that committing it
   gives this commit a real `push`-triggered CI run on its own exact SHA --
   closing the gap `check_pretag_ci.py` correctly refused to let slide.

## Resulting chain

- M (release PR #73 merge): `4514cfc` -- own CI green (run `36880132412` on
  the PR head `a5c5463`; `pr-policy-gate` fires post-merge).
- T (sweep evidence): `bb91f5a` -- own CI green (run `36881903928`); sweep
  `release`+`readme` SUCCESS, e2e 23/23 cells.
- Board row: `7718c3a` -- `pathway_forward.md` Tasks-done row for v5.21.0,
  inherited into every descendant's tree.
- This commit -- carries the row (inherited from `7718c3a`) and this log
  (non-ignored), so it is the first SHA in this chain able to satisfy
  `check_pretag_ci.py` end to end: own CI green, own `publish.yml` green,
  sufficient sweep evidence, and release docs ready, all on the SAME sha.

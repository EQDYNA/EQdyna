# History rewrite, 2026-09-30

Owner-approved one-time exception to the no-history-rewrite rule (board item 36): `scec_archive/` was committed and removed on 2026-09-16 but stayed in history, 233 MB of a 305 MB clone. It was dropped from every commit with `git filter-repo --path scec_archive/ --invert-paths`.

- **Result:** a fresh clone went from ~305 MB to 74 MB (`.git`), 20 s.
- **Unchanged:** every file tree. All 101 branch and tag tips have the same tree hash before and after; no tip contained `scec_archive/`. The archive itself is untracked and was never touched (894/894 checksums OK after).
- **Changed:** every commit SHA from 2026-09-16 on. All 42 tags were force-updated to their rewritten commits; GitHub Releases follow the tag names.
- **Old SHAs** cited anywhere in docs, the board, session logs or commit messages written before this date refer to the old history. Translate them with `COMMIT_MAP_2026-09-30.txt` (`old new`, full SHAs, 1168 commits): `grep ^<old-prefix> docs/notes/COMMIT_MAP_2026-09-30.txt`.
- **Evidence:** the seven `docs/evidence/sweep-<sha>/summary.json` files were renamed to their new short SHA; `sha` holds the new SHA and `sha_pre_rewrite` the old one.
- **Backup:** a full pre-rewrite mirror and a verified bundle (all refs, including GitHub PR refs) are kept off-repo by the owner.
- **Existing clones** must be re-cloned (or `git fetch --force` + `git reset --hard origin/master` on a checkout with no local work).

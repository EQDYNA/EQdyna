# scripts/scec — re-fetch EQdyna's own SCEC/USGS submissions

This is the TOOLING. The data it fetches is not tracked and must not be.

## Why this is here and not in scratch/

`pathway_forward.md` item 28: TPV29's cross-code validation numbers (median
rupture-time differences 0.076 / 0.028 / 0.006 s at 500 / 200 / 100 m vs
EQdyna's 2015 submission, ruptured area 749-756 km², Mw 7.45) were prose-only,
and the baseline they were measured against lived in gitignored `scratch/`
along with the tooling that produced it. Nobody could reproduce them from a
clean clone — a rule 4 violation, since a number has to travel with something
that regenerates it.

The data is **484 MB** and correctly stays out of git. The tooling is **116 KB**
and had no business being untracked. That is the whole fix: track what
regenerates the baseline, not the baseline.

Pure standard library — no new CI dependency (checked against
`.github/workflows/test.yml`'s pip line, and `test_ci_dependencies.py` now
scans this directory).

## Layout

    fetch_all_dliu.py   fetch every dliu* submission from the cvws portal
    organize.py         ARCHIVE/PAGES/CGI paths, RMS matrix parsing
    make_index.py       regenerate scec_archive/INDEX.md
    provenance.py       per-submission provenance stamps
    *.json              enumeration / manifest / layout / index of what the
                        portal is known to hold

`organize.py` resolves `ARCHIVE` as `<this dir>/../../scec_archive`, which is
the repo root either from here or from the original `scratch/scec_archive/` —
the move preserved the depth deliberately.

## What is archived, and what it is for

Published EQdyna submissions, fetched 2026-09-14 from
`https://strike.scec.org/cvws/`: tpv29, tpv30, tpv34, tpv35, tpv36, tpv37,
tpv104, tpv105-3d. These are the ONLY independent validation targets this repo
has — every `test.reference.results/` artifact is EQdyna compared against
itself, so it catches regression but cannot catch a shared error.

Published standing is worth reading before trusting a coarse run:

| benchmark | resolution | median RMS rupture time | rank |
|---|---|---|---|
| tpv29 | 100 m | 87.2 ms | 15/15 |
| tpv29 | 50 m | 33.0 ms | 12/15 |
| tpv30 | 100 m | 287.7 ms | 14/14 |
| tpv30 | 50 m | 120.6 ms | 12/14 |
| tpv30 | 25 m | 55.1 ms | 8/14 |
| tpv34 | 50 m | 32.4 ms | 13/14 |

The trend is the point: coarse runs ranked LAST and improved monotonically
with refinement. The 500 m spacings the e2e sweep gates at are regression
checks, not physics claims, and nothing in this repo should be read as a
statement about accuracy at that resolution.

## Re-fetching

    python3 scripts/scec/fetch_all_dliu.py       # writes scec_archive/
    python3 scripts/scec/make_index.py           # regenerates INDEX.md

Depends on the public portal being reachable. If it moves or retires, this
tooling is the record of what was there and when.

## Verifying the on-disk archive against a tracked manifest

`CHECKSUMS.sha256` (693 lines, ~80 KB, standard `sha256sum` format, paths
relative to `scec_archive/`) is tracked so the 484 MB archive's integrity is
auditable without putting the archive itself in git — "it is an asset, hash
but no need to be in git" (owner, 2026-09-16, after `f2c9851` accidentally
committed all 484 MB: a `.gitignore` inline comment on the `scec_archive/`
line made the whole line an unmatchable literal pattern, since gitignore has
no inline-comment syntax; `git rm --cached -r scec_archive/` untracked it
going forward, `git gc` reclaimed the loose objects, and this manifest is the
integrity check that replaces "it's in git so I'd notice if it changed").

    python3 scripts/scec/organize.py --verify

Reports OK / MISMATCH / MISSING / EXTRA counts and exits non-zero on any of
the three — the manifest is meant to enumerate the full archive, so an EXTRA
file (present on disk, absent from the manifest, e.g. after a re-fetch that
grew the archive before this file was regenerated) is drift worth seeing,
not something to pass silently (rule 2). Regenerate the manifest deliberately
when the archive legitimately changes; do not loosen this check instead.

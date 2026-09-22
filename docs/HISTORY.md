# EQdyna legacy history, 2014-2016

The first 8 commits of this repository are version drops of code written
between **2014-09-12 and 2016-10-05**, bulk-imported on **2018-11-26**. Git
stamps all 8 with the import date, so `git log`, `git blame`, GitHub's
history view and any date-ordered analysis see a flat wall on one day. The
real dates were never lost -- they are written in the commit MESSAGES -- but
no tool reads them.

This file is that timeline, parsed out of those messages.

**It is a record, not a repair.** Restoring the true author dates would
rewrite every downstream SHA (724 commits, 33 tags, both forks);
item 36 decided against a history rewrite and that stands. `git notes` was
the other cheap option and was rejected: notes live in a ref most clones
never fetch and most viewers never show.

**Nothing below is inferred.** Every cell is parsed from the commit it
describes by `testsys/regression/test_history_table.py`, which also VERIFIES
this table on every regression run -- if a message and this table ever
disagree, the test fails. Two commits state only a month, so the precision
column says `month` rather than inventing a day. One declares no version at
all (`221b3ad` mentions 3.2.1 only as the version its MPI work was merged
INTO), so its version reads `(unstated)`.

One detail the earlier prose account got slightly wrong, corrected here from
the log: the eight share an AUTHOR date of 2018-11-26, but `d1dc019`'s
COMMITTER date is 2018-11-27 -- the import ran past midnight. This table keys
on the author date, and shows both stamps.

Regenerate with `python3 testsys/regression/test_history_table.py --regenerate`.

## The pre-2018 source itself

A legacy archive of the pre-2018 source exists on disk at
`scratch/legacy_import/Examples` (gitignored like every `scratch/` asset),
hash-manifested at `docs/evidence/legacy_import_archive.sha256`. It is NOT
needed to read this table: the dates come from the commit messages, which are
in the repo. The archive's own bundled `.git` covers only 4.1.1/4.1.2/4.1.3
and stamps all three 2018-11-26 -- a second snapshot import, not a timeline,
and importing it would add a second flat wall rather than fix the first.

## The timeline

<!-- BEGIN GENERATED TABLE -- testsys/regression/test_history_table.py -->

| # | commit | TRUE authored date | precision | author(s) | version declared | git stamp (author / committer) | the commit's own verification line (excerpt) |
|---|--------|--------------------|-----------|-----------|------------------|-------------------------------|----------------------------------------------|
| 1 | `d8d65d0` | **2014-09-12** | day | Benchun Duan | 3.1 | 2018-11-26 / 2018-11-26 | -- |
| 2 | `0a983b7` | **2015-09** | month | Dunyu Liu | 3.1.1 | 2018-11-26 / 2018-11-26 | verified against SCEC TPV29&30 |
| 3 | `16d334c` | **2015-09** | month | Dunyu Liu | 3.1.2 | 2018-11-26 / 2018-11-26 | -- |
| 4 | `7d603a1` | **2015-09-19** | day | Dunyu Liu | 3.2.1 | 2018-11-26 / 2018-11-26 | verified against the model with PML and Q model in Ma |
| 5 | `221b3ad` | **2016-09-29** | day | Dunyu Liu, Bin Luo | (unstated) | 2018-11-26 / 2018-11-26 | verified against SCEC TPV8 |
| 6 | `3d512a8` | **2016-10-04** | day | Dunyu Liu | 4.1.1 | 2018-11-26 / 2018-11-26 | verified against SCEC TPV8 |
| 7 | `f2c7c0e` | **2016-10-04** | day | Dunyu Liu | 4.1.2 | 2018-11-26 / 2018-11-26 | verified against SCEC TPV8 |
| 8 | `d1dc019` | **2016-10-05** | day | Dunyu Liu, Bin Luo | 4.2 | 2018-11-26 / 2018-11-27 | verified against SCEC TPV104 |

<!-- END GENERATED TABLE -->

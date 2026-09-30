#! /usr/bin/env python3
"""ONE copy of the two release-readiness text checks, operating on TEXT
(never a path) so both the POST-hoc reader
(`testsys/regression/test_release_complete.py`'s `check_pathway_tasks_done_row`
/ `check_readme_news_leads_with_this_version`, which read the live working
tree) and the PRE-hoc gate (`testsys/regression/check_pretag_ci.py`'s
`evaluate_release_docs_ready`, which reads a candidate commit's OWN TREE via
`git show <sha>:<path>`, never the live checkout) share the identical regex
and verdict logic. Two copies of this drifting apart -- one file's idea of
"has a Tasks-done row" disagreeing with the other's -- would silently
reopen exactly the ordering gap rule 15 step 4 / the owner's (b)/(d)
requirement exist to close.

WHY THE POST-HOC AND PRE-HOC HALVES BOTH NEED THIS. The owner's requirement
(b): "everything a post-tag check reads (board Tasks-done row, notes,
Release draft) lands before the tag." `check_pretag_ci.py` is what enforces
that BEFORE `git tag` runs; `test_release_complete.py` is what re-confirms
it AFTER, against the tagged tree, same as it already does for CI-green and
sweep-evidence. A regex change that only reached one side would make the
pre-tag gate and the post-hoc regression check disagree about the same
release -- passing one, failing the other -- which is a worse failure mode
than either check being wrong by itself.
"""
import re

TASKS_DONE_ROW_MENTION_RE = r'v%s\b'
TASKS_DONE_ROW_DATED_RE = r'^\| 20\d\d-\d\d-\d\d \|.*v%s'
README_NEWS_LEAD_RE = re.compile(r'^\* \d{8} v([0-9.]+) release notes', re.M)


def tasks_done_row_present(pathway_text, v):
    """(ok, message) -- does `pathway_text` (pathway_forward.md's content,
    from ANY tree: live working copy or `git show sha:pathway_forward.md`)
    carry a dated Tasks-done row naming version `v` (e.g. '5.20.2', no
    leading 'v')? Two sub-checks, same as the original
    check_pathway_tasks_done_row: a bare mention is not enough, it must be
    a '| YYYY-MM-DD | ... vX.Y.Z' row (rule 15 step 4)."""
    if not re.search(TASKS_DONE_ROW_MENTION_RE % re.escape(v), pathway_text):
        return False, 'never mentions v%s' % v
    if not re.search(TASKS_DONE_ROW_DATED_RE % re.escape(v), pathway_text, re.M):
        return False, ('mentions v%s but has no dated Tasks-done row for it '
                       '(rule 15 step 4)' % v)
    return True, 'has a dated Tasks-done row for v%s' % v


def readme_news_leads_with(readme_text, v):
    """(ok, message) -- does `readme_text`'s LEADING '* YYYYMMDD vX.Y.Z
    release notes' block already name version `v`? Same regex and rule
    (rule 15 step 3: the current release leads README, the previous one
    moves to pastReleaseNotes.md) as the original
    check_readme_news_leads_with_this_version."""
    m = README_NEWS_LEAD_RE.search(readme_text)
    if not m:
        return False, ("has no '* YYYYMMDD vX.Y.Z release notes' block "
                       "(rule 15 step 3)")
    if m.group(1) != v:
        return False, ('leading News block is v%s, not v%s (rule 15 step 3: '
                       'the current release leads README)' % (m.group(1), v))
    return True, 'leads with v%s' % v

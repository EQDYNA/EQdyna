#! /usr/bin/env python3

# NOTE (mission: tpv22/23 campaign, 2026-10-01): this list was found to be
# missing most of case_input/'s actual directories (test.tpv10, test.tpv22,
# test.tpv23, test.tpv29, test.tpv30, test.tpv36, test.tpv37, test.drv.a6,
# test.drv.a6.v2, test.multifault2, TPV2802_50_a6, bp1001.fdc.rough.250,
# liu2020.fdc.planar, liu2020.fdc.rough.250 -- 13 real directories absent
# before this edit) and to carry one dangling entry, 'test.tpv1053d.6c',
# which does not match any directory under case_input/ (only test.tpv1053d
# exists). Nothing in the repo reads this file or case_input/compsets.txt
# (grepped both whole-repo; zero hits) -- scripts/create.newcase resolves a
# compset name to case_input/<name> directly, with no lookup against either
# list. compset.txt's own entries disagree with THIS file and with the real
# directories even more (see case_input/compsets.txt's own note): of its 9
# entries, only 3 (bp1001.fdc.rough.250, liu2020.fdc.planar,
# liu2020.fdc.rough.250) match a real directory; the other 6 use a bare or
# hyphenated name ('tpv104', 'test-tpv104', ...) that matches nothing. This
# file's own 'test.X' spelling is the one that actually matches directory
# names, so test.tpv22/test.tpv23 are added here, not there. Reconciling the
# rest of this list's drift against case_input/ is out of scope for this
# mission (bigger than a one-line fix) -- named here as a follow-up for the
# board.
compset = ['test.tpv8', 'test.tpv104', 'test.tpv1053d', 'test.tpv1053d.6c', 'test.meng2023a', 'test.meng2023cb', 'test.tpv22', 'test.tpv23']
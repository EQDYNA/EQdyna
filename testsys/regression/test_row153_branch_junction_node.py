#! /usr/bin/env python3
"""
Regression guard (row 153 checkpoint 2, physics defect): test.tpv24's branch
fault never ruptured (measured 0/160 deduped fault-2 nodes by 5s, 1/170 even
at the case's own 12s term), while Barall's independent FaultMod 100m
reference (SCEC CVWS tpv24 submission) ruptures the branch fully, first at
t=3.08s.

Root cause, found by measurement (falsifying three other suspects first --
see the fix commit message for the full account): the near-junction branch
node (BRANCH_FXMIN, 1 dx = 1 km from the junction) was locked
mu_s=UNBREAKABLE_FS=1000 in build_on_fault_vars_branch, on the theory that
it was "a true border of fault 2... slip goes to zero at the junction, like
any other border" (spec p.5). It is not a border: the junction point itself
(x0 = fxmin - dx = 0) is EXCLUDED from fault 2's box entirely (owned by
fault 1 only). The near-junction node is the first node *inside* the
branch -- and per Barall's branchst010dp100 station (L=1 km from the
junction), it is the node that ruptures FIRST, exactly where the main
fault's rupture is supposed to jump onto the branch. Locking it unbreakable
made that structurally impossible regardless of any correct stress/frame
setup underneath it.

This check is deliberately NOT a stress-value check (the initial branch
shear/normal stress was already correct before this fix -- verified against
Barall's branchst*dp* t=0 columns, see the fix commit). A stress-value
assertion alone would NOT have caught this regression. What catches it is
asserting the near-junction node's static friction coefficient is the
ordinary breakable MU_S, not the unbreakable border value -- and, as a
second line of defense, checking the resolved initial branch traction
against Barall's own measured t=0 reference values, so a future change to
the stress derivation is caught independently of the friction-field check.

MUTATION TEST (both-ways, rule 14a): reverting the fix (restoring
`on_left_edge = (ix == 0)` and `if on_left_edge or on_right_edge or
on_bottom_edge:`) makes check_junction_node_breakable fail (mu_s reads
UNBREAKABLE_FS=1000.0, not MU_S=0.18) -- demonstrated by hand against the
pre-fix commit, this session.

No solver run needed -- this imports tpv24_25_common directly and inspects
its own on_fault_vars arrays, so it runs in well under a second and needs
no built bin/eqdyna.
"""
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CASE_INPUT = os.path.join(REPO_ROOT, 'case_input', 'test.tpv24')
SCRIPTS = os.path.join(REPO_ROOT, 'scripts')

# Barall FaultMod 100m reference, branchst010dp100 (L=1 km along branch
# strike from the junction, 10 km down-dip), t=0 row:
#   h-shear-stress =  2.420747E+01 MPa, n-stress = -1.456933E+02 MPa
# (scratch dataset path deliberately not hard-coded here -- read-only,
# external, not committed; the numeric reference below is transcribed from
# it, same convention as every other spec-derived constant in this repo).
BARALL_JUNCTION_H_SHEAR_MPA = 24.20747
BARALL_JUNCTION_N_STRESS_MPA = -145.6933
TOL_MPA = 0.01  # ~5 significant figures, matches the gate's own dx=1000m coarseness


def _load_common():
    sys.path.insert(0, SCRIPTS)      # defaultParameters, lib (tpv24_25_common's own deps)
    sys.path.insert(0, CASE_INPUT)
    import tpv24_25_common as common  # noqa
    return common


def check_junction_node_breakable():
    common = _load_common()
    nfz = round((common.FZMAX - common.FZMIN) / 1000.0 + 1)
    fz = [common.FZMIN + i * 1000.0 for i in range(nfz)]
    nfx2 = round((common.BRANCH_FXMAX - common.BRANCH_FXMIN) / 1000.0 + 1)
    fx2 = [common.BRANCH_FXMIN + i * 1000.0 for i in range(nfx2)]
    # TPV24 coefficients (user_defined_params.py's own values)
    v = common.build_on_fault_vars_branch(fx2, fz, nfx2, nfz,
                                           b22=0.926793, b33=1.073206, b23=-0.169029)
    iz10km = fz.index(-10000.0)
    ix_junction = 0  # BRANCH_FXMIN, one dx from the excluded junction column
    mu_s = v[iz10km, ix_junction, 1]
    assert mu_s == common.MU_S, (
        'near-junction branch node (x=%.0f, z=-10000) has mu_s=%r, expected the '
        'ordinary breakable MU_S=%r -- it has been re-locked unbreakable '
        '(UNBREAKABLE_FS=%r), which froze the branch before (row153 ckpt2 '
        'regression)' % (common.BRANCH_FXMIN, mu_s, common.MU_S, common.UNBREAKABLE_FS))

    shear_mpa = v[iz10km, ix_junction, 8] / 1.0e6
    norm_mpa = v[iz10km, ix_junction, 7] / 1.0e6
    assert abs(shear_mpa - BARALL_JUNCTION_H_SHEAR_MPA) < TOL_MPA, (
        'near-junction branch initial shear stress %.5f MPa does not match Barall '
        'reference %.5f MPa (branchst010dp100, t=0)' % (shear_mpa, BARALL_JUNCTION_H_SHEAR_MPA))
    assert abs(norm_mpa - BARALL_JUNCTION_N_STRESS_MPA) < TOL_MPA, (
        'near-junction branch initial normal stress %.5f MPa does not match Barall '
        'reference %.5f MPa (branchst010dp100, t=0)' % (norm_mpa, BARALL_JUNCTION_N_STRESS_MPA))

    # far (true) edge and bottom must STILL be locked -- the fix must not
    # have swung the other way and unlocked every border.
    far_edge_mu_s = v[iz10km, nfx2 - 1, 1]
    bottom_mu_s = v[0, ix_junction, 1]
    assert far_edge_mu_s == common.UNBREAKABLE_FS, (
        'far branch edge (x=%.0f) must stay unbreakable (true free edge), got mu_s=%r'
        % (common.BRANCH_FXMAX, far_edge_mu_s))
    assert bottom_mu_s == common.UNBREAKABLE_FS, (
        'bottom branch edge (z=%.0f) must stay unbreakable (true free edge), got mu_s=%r'
        % (common.FZMIN, bottom_mu_s))


def main():
    print('Regression guard: row 153 ckpt2 -- branch junction node must stay breakable')
    rc = 0
    try:
        check_junction_node_breakable()
    except AssertionError as e:
        print('  FAIL  check_junction_node_breakable: %s' % e)
        rc = 1
    else:
        print('  PASS  check_junction_node_breakable')
    print('\n%s test_row153_branch_junction_node' % ('FAIL' if rc else 'SUCCESS'))
    return rc


if __name__ == '__main__':
    sys.exit(main())

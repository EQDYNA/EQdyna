#! /usr/bin/env python3
"""
Regression guard: every on-fault station of test.tpv22/test.tpv23 sits on a
node of ITS fault (rules 2, 10).

THE DEFECT: SCEC TPV22/23 Part 6 p.11 lists fault #2 station 8 at 6.5 km along
strike. Fault #2's nodes are at -5000 + k*dx, so at TPV22's dx=200 m 6.5 km is
midway between two nodes; setOnFaultStation (meshgen.f90) and its port match x
within tol=1e-5 m, so the station matched nothing and no
faultstft2_*dp100 file was written for it. The spec (p.12) permits moving a
station to the nearest node; tpv22_23_common.snap_to_node does that.

Pure arithmetic on build_params(); no solver, no MPI. Exits non-zero on failure.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path[:0] = [os.path.join(ROOT, 'scripts'), os.path.join(ROOT, 'case_input', 'test.tpv22')]
TOL_M = 1.0e-5  # globalvar.f90 module tol, as meshgen.f90/meshgen.py use

import tpv22_23_common as c  # noqa: E402


def main():
    fail = []
    n = 0
    for tpv, z2 in ((22, -1600.0), (23, 1000.0)):
        par = c.build_params(tpv, z2)
        for ift, (xr, sts) in enumerate(zip((c.FAULT1_XRANGE, c.FAULT2_XRANGE),
                                            par.st_coor_on_fault_per_fault), 1):
            for x_km, z_km in sts:
                n += 1
                for name, v, lo, d in (('x', x_km * 1000.0, xr[0], par.dx),
                                       ('z', z_km * 1000.0, c.FZMIN, par.dz)):
                    off = abs((v - lo) / d - round((v - lo) / d)) * d
                    if off > TOL_M:
                        fail.append('tpv%d fault %d station (%g,%g) km: %s is %.3f m '
                                    'off the node grid' % (tpv, ift, x_km, z_km, name, off))
    # the spec's own 6.5 km must have moved to its nearest node at dx=200 and
    # stayed put at dx=250 (both directions of the guard)
    s22 = c.build_params(22, -1600.0).st_coor_on_fault_per_fault[1][7]
    s23 = c.build_params(23, 1000.0).st_coor_on_fault_per_fault[1][7]
    if s22 != [6.6, -10.0] or s23 != [6.5, -10.0]:
        fail.append('station 8 expected [6.6,-10] (tpv22) / [6.5,-10] (tpv23), got %r / %r' % (s22, s23))
    if fail:
        print('FAIL\n' + '\n'.join(fail))
        return 1
    print('SUCCESS: %d on-fault stations checked (tpv22+tpv23, both faults), all on-grid' % n)
    return 0


if __name__ == '__main__':
    sys.exit(main())

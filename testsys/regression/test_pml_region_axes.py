#! /usr/bin/env python3
"""
Regression guard (PROJECT_RULES.md rule 10) for the PML region-14 bug
fixed 2026-09-13: the region cascade compared the y coordinate against
the x-axis bound (y >= xmax0 / xc(2) > xmax2) in all four copies of the
cascade, corrupting damping in the -x/+y PML corner (moved tstk by up
to 29 MPa in test.drv.a6). This guard fails if any cross-axis
coordinate/bound comparison reappears in the PML sources.
"""
import os, re, sys

root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
bad = []
# y-coordinate (y or xc(2)) vs x bound; x-coordinate vs y bound; z vs x/y bound
patterns = [r'\by\s*[<>]=?\s*x(?:max|min)', r'\bxc\(2\)\s*[<>]=?\s*x(?:max|min)',
            r'\bx2\s*[<>]=?\s*y(?:max|min)', r'\bxc\(1\)\s*[<>]=?\s*y(?:max|min)',
            r'\bxc\(3\)\s*[<>]=?\s*[xy](?:max|min)', r'\bz\s*[<>]=?\s*[xy](?:max|min)']
for fname in ['src/comdampv.f90', 'src/assembleGlobalKU.f90']:
    for n, line in enumerate(open(os.path.join(root, fname)), 1):
        code = line.split('!')[0]
        for p in patterns:
            if re.search(p, code):
                bad.append(f'{fname}:{n}: {line.strip()}')
if bad:
    print('FAIL test_pml_region_axes: cross-axis PML comparison found')
    for b in bad:
        print('  -', b)
    sys.exit(1)
print('SUCCESS test_pml_region_axes')
